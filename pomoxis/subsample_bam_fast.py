import argparse
import functools
import logging
import multiprocessing import cpu_count, Process, Queue, Array
import os

from intervaltree import IntervalTree, Interval
import numpy as np
import pysam


from pomoxis.util import parse_regions, Region
from pomoxis.coverage_from_bam import coverage_summary_of_region
from pomoxis.stats_from_bam import stats_from_aligned_read


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser(
        prog='subsample_bam',
        description='Quickly subsample a bam to uniform depth',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam',
        help='input bam file.')
    parser.add_argument('depth', nargs='+', type=int,
        help='Target depth.')
    parser.add_argument('-o', '--output_prefix', default='sub_sampled',
        help='Output prefix')
    parser.add_argument('-r', '--regions', nargs='+',
        help='Only process given regions.')
    parser.add_argument('-p', '--profile', type=int, default=1000,
        help='Stride in genomic coordinates for depth profile.')
    parser.add_argument('-O', '--orientation', choices=['fwd', 'rev'],
        help='Sample only forward or reverse reads.')
    parser.add_argument('-t', '--threads', type=int, default=-1,
        help='Number of threads to use.')
    parser.add_argument('-q', '--quality', type=float,
        help='Filter reads by mean qscore.')
    parser.add_argument('-a', '--accuracy', type=float,
        help='Filter reads by accuracy.')
    parser.add_argument('-c', '--coverage', type=float,
        help='Filter reads by coverage (what fraction of the read aligns).')
    parser.add_argument('-l', '--length', type=int, default=None,
        help='Filter reads by read length.')

    parser.add_argument('--output-fasta', default=None,
        help='Create a FASTA file of the subsampled reads.')

    eparser = parser.add_mutually_exclusive_group()
    eparser.add_argument('--any_fail', action='store_true',
        help='Exit with an error if any region has insufficient coverage.')
    eparser.add_argument('--all_fail', action='store_true',
        help='Exit with an error if all regions have insufficient coverage.')

    uparser = parser.add_argument_group('Uniform sampling options')
    uparser.add_argument('-x', '--patience', default=5, type=int,
        help='Maximum iterations with no change in median coverage before aborting.')
    uparser.add_argument('-s', '--stride', type=int, default=1000,
        help='Stride in genomic coordinates when searching for new reads. Smaller can lead to more compact pileup.')

    args = parser.parse_args()
    if args.threads == -1:
        args.threads = cpu_count()

    with pysam.AlignmentFile(args.bam) as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))

        if args.regions is not None:
            regions = parse_regions(args.regions, ref_lengths=ref_lengths)
        else:
            regions = [Region(ref_name=r, start=0, end=ref_lengths[r]) for r in bam.references]

    enough_depth = []

    processes = []
    work_queue = Queue()
    ret_read_q = Queue()
    ret_bam_q = Queue()

    for _ in range(args.threads):
        p = Process(target=subsample_region_uniformly, args=(work_queue, ret_bam_q, ret_read_q, args))
        process.append(p)
    for p in processes:
        p.start()

    if args.output_fasta:
        p = Process(target=write_reads, args=(ret_read_q, args))
        processes.append(p)

    for region in regions:
        work_queue.put({
            "region": region,
        })

    # Add sentinels to queue
    for _ in range(args.threads):
        work_queue.put(None)

    # Wait for threads to finish processing regions
    for p in processes:
        p.join()


    # Write the BAM out (it's mp-ready, but uses one thread for now)
    post_processes = []
    p = Process(target=write_bam, args=(ret_bam_q, args))
    post_processes.append(p)
    ret_bam_q.put(None)
    p.join()



    if args.any_fail and not all(enough_depth):
        raise RuntimeError('Insufficient read coverage for one or more requested regions.')
    if args.all_fail and all([not x for x in enough_depth]):
        raise RuntimeError('Insufficient read coverage for all requested regions.')


def filter_read(r, bam, args, logger):
    """Decide whether a read should be filtered out, returning a bool"""

    # primary alignments
    if (r.is_secondary or r.is_supplementary):
        return True

    # filter orientation
    if (r.is_reverse and args.orientation == 'fwd') or \
        (not r.is_reverse and args.orientation == 'rev'):
        return True

    # filter quality
    if args.quality is not None:
        mean_q = np.mean(r.query_qualities)
        if mean_q < args.quality:
            logger.debug("Filtering {} by quality ({:.2f}).".format(r.query_name, mean_q))
            return True

    # filter accuracy or alignment coverage
    if args.accuracy is not None or args.coverage is not None or args.length is not None:
        stats = stats_from_aligned_read(r, bam.references, bam.lengths)
        if args.accuracy is not None and stats['acc'] < args.accuracy:
            logger.info("Filtering {} by accuracy ({:.2f}).".format(r.query_name, stats['acc']))
            return True
        if args.coverage is not None and stats['coverage'] < args.coverage:
            logger.info("Filtering {} by coverage ({:.2f}).".format(r.query_name, stats['coverage']))
            return True
        if args.length is not None and stats['read_length'] < args.length:
            logger.info("Filtering {} by length ({:.2f}).".format(r.query_name, stats['length']))
            return True
    # don't filter
    return False


def get_cursor_reads(bam, cursor, max_abs_dist=500, seen=[]):
    # Pileup reads around the cursor and organise them into a structure
    # based on their start distance compared to the cursor
    candidates_by_distance = {}

    # Call pileup at cursor
    #   * cursor is 0 indexed
    #   * stepper "all" will drop secondary reads etc.
    #   * min_base_quality is 0, if you want to drop reads based on quality, let filter_read do it
    for pcol in bam.pileup(contig=region.ref_name, start=cursor, end=cursor+1, stepper="all", min_base_quality=0):
        for read in pcol.pileups:
            if filter_read(read, args):
                continue

            # drop read if was seen recently
            if read.query_name in seen:
                continue

            # Select reads that begin before, or ON the cursor
            if read.reference_start == cursor+1:
                break

            delta = abs(read.reference_start - cursor) # force positive as we're only looking left of cursor now

            # limit max distance from cursor
            if abs(delta) > max_abs_dist:
                continue

            if delta not in candidates_by_distance:
                candidates_by_distance[delta] = []
            candidates_by_distance[delta].append(read)

    return candidates_by_distance

def select_cursor_reads(d_reads, n=1):
    selected = []
    bins = sorted(d_reads.keys())
    bins_round = 0

    while len(selected) < n:
        # Keep adding bins until at least N bins have been observed
        selected.extend( d_reads[bins[bins_round]] )

        bins_round += 1
        if len(selected) < n and bins_round == len(bins):
            logger.warning("[select_cursor_reads] Eligible reads underfilled. Desired %d, Acquired %d" % (n, len(selected)))
            return selected # return what we have without randomising

    return np.random.choice(selected, n, replace=False) # select N random reads without replacement


def write_reads(read_q, args):
    reads_fh = open(args.output_fasta, 'w')
    dead_threads = 0

    while True:
        read = read_q.get()
        if read is None:
            dead_threads += 1
            if dead_threads == args.threads:
                break
            continue

    reads_fh.close()



def subsample_region_uniformly(work_q, bam_ret_q, read_ret_q, args):
    while True:
        work = work_q.get()
        if work is None:
            # If sentinel, push sentinel to return queue and finish
            read_ret_q.put(None)
            return
        else:
            region = work["region"]

        logger = logging.getLogger(region.ref_name)
        logger.info("Pulled %s from queue" )

        # Begin cursor at region start
        cursor = region.start

        # Initialise a coverage track for each layer of desired depth
        track_ends = np.full((1, args.depth[0]), cursor)

        # Keep track of reads and coverage
        n_reads = 0
        seen_reads = set([])
        track_reads = {}
        track_cov = np.zeros(region.end - region.start)

        while cursor < region.end:
            # Count the number of tracks that need to be filled at this position
            tracks_at_cursor = np.argwhere(track_ends == cursor).flatten()

            # Gather the reads, then randomly select one for each track
            eligible_reads = get_cursor_reads(bam, cursor, seen=seen_reads)
            chosen_reads = select_cursor_reads(eligible_reads, len(tracks_at_cursor))

            for read_i, read in enumerate(chosen_reads):
                n_reads += 1
                curr_track = tracks_at_cursor[read_i]
                track_reads[curr_track] = read
                track_ends[curr_track] = read.reference_end
                track_cov[read.reference_start - region.start : read.reference_end - region.start] += 1

                bam_ret_q.put(read.query_name) # send read name to write to new BAM after all threads are done
                seen_reads.add(read.query_name)

                if args.output_fasta:
                    read_ret_q.put(read.to_string()) # send read data to be written to fasta/fastq

            # Move cursor to earliest track end
            curr_track = np.argmin(track_ends)
            cursor = track_ends[curr_track]

        median_depth = np.median(track_cov)
        stdv_depth = np.std(track_cov)
        logger.info(u'region: {}, reads: {}, target: {}, depth: {:.0f}X (\u00B1{:.1f}).'.format(region.ref_name, n_reads, args.depth[0], median_depth, stdv_depth))


def _nearest_overlapping_point(src, point):
    """Find the interval with the closest start point to a given point.

    :param src: IntervalTree instance.
    :param point: query point.

    :returns: Interval instance of interval with closest start.

    """
    items = src.at(point)
    if len(items) == 0:
        return None
    items = sorted(items, key=lambda x: x.end - x.begin, reverse=True)
    items.sort(key=lambda x: abs(x.begin - point))
    return items[0]


def _write_bam(bam, prefix, region, sequences):
    # filtered bam
    sequences = set(sequences)
    taken = set()
    output = '{}_{}.{}'.format(prefix, region.ref_name, os.path.basename(bam))
    src_bam = pysam.AlignmentFile(bam, "rb")
    out_bam = pysam.AlignmentFile(output, "wb", template=src_bam)
    for read in src_bam.fetch(region.ref_name, region.start, region.end):
        if read.query_name in sequences and read.query_name not in taken:
            out_bam.write(read)
            taken.add(read.query_name)
    src_bam.close()
    out_bam.close()
    pysam.index(output)


def _write_coverage(prefix, region, coverage, profile):
    # depth profile
    output = '{}_{}.depth'.format(prefix, region.ref_name)
    end = profile * (len(coverage) // profile)
    cov_blocks = coverage[0:end].reshape(-1, profile)
    depth_profile = np.mean(cov_blocks, axis=1, dtype=np.uint32)
    start = region.start + profile // 2
    positions = (start + profile * x for x in range(len(depth_profile)))
    with open(output, 'w') as fh:
        fh.write("position\tdepth\n")
        for pos, depth in zip(positions, depth_profile):
            fh.write("{}\t{}\n".format(pos, depth))


if __name__ == '__main__':
    main()
