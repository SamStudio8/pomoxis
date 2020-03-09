import argparse
import functools
import logging
import os
from multiprocessing import cpu_count, Process, Queue, Manager

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
    parser.add_argument('depth', type=int,
        help='Target depth.')
    parser.add_argument('-o', '--output', required=True, default='-',
        help='path to output subsampled BAM [default -]')
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
    read_ret_q = Queue()

    manager = Manager()
    bam_ret_d = manager.dict()

    for _ in range(args.threads):
        p = Process(target=subsample_region_uniformly, args=(work_queue, bam_ret_d, read_ret_q, args))
        processes.append(p)
    if args.output_fasta:
        p = Process(target=write_reads, args=(read_ret_q, args))
        processes.append(p)
    pp = Process(target=write_bam, args=(bam_ret_d, args))
    processes.append(pp)
    for p in processes:
        p.start()

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
    pp.join()


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

        # Write read (by parsing pysam AlignedSegment as str)
        fields = read.split()
        reads_fh.write(">%s\n%s\n" % (fields[0], fields[9]))

    reads_fh.close()


def write_bam(uuid_d, args):
    src_bam = pysam.AlignmentFile(args.bam, "rb")
    out_bam = pysam.AlignmentFile(args.output, "wb", template=src_bam)

    # One-pass through BAM means the output should be sorted?
    for read in src_bam.fetch():
        if "%s:%d" % (read.query_name, read.reference_start) in uuid_d:
            out_bam.write(read)
    src_bam.close()
    out_bam.close()


def subsample_region_uniformly(work_q, bam_ret_d, read_ret_q, args):
    while True:
        work = work_q.get()
        if work is None:
            # If sentinel, push sentinel to return queue and finish
            read_ret_q.put(None)
            return
        else:
            region = work["region"]

        logger = logging.getLogger(region.ref_name)
        logger.info("Pulled %s from queue" % region.ref_name)

        bam = pysam.AlignmentFile(args.bam)

        # Begin cursor at region start
        CURSOR_DIST = 500
        closest_cursor = region.start
        if closest_cursor == 0:
            closest_cursor += CURSOR_DIST
        cursors = {closest_cursor: []}

        # Initialise a coverage track for each layer of desired depth
        track_ends = np.full(args.depth, max(closest_cursor, CURSOR_DIST))

        # Keep track of reads and coverage
        n_reads = 0
        seen_reads = set([])
        track_reads = {}
        track_cov = np.zeros(region.end - region.start)



        for read_i, read in enumerate(bam.fetch(contig=region.ref_name, start=region.start, end=region.end)):
            if filter_read(read, bam, args, logger):
                continue

            # drop read if used already
            if read.query_name in seen_reads:
                continue

            # If in range of at least one cursor, randomly pick a close cursor to give this read to
            if read.reference_start >= closest_cursor-CURSOR_DIST and read.reference_start < closest_cursor:
                eligible_keys = []
                keys = list(cursors.keys())
                for c_i in np.argsort(keys):
                    c_cursor = keys[c_i]
                    if read.reference_start >= c_cursor-CURSOR_DIST and read.reference_start < c_cursor:
                        eligible_keys.append(c_i)
                #print(cursors.keys(), closest_cursor, read.reference_start)
                #print({key: len(cursors[key]) for key in cursors})
                rand_i = np.random.choice(eligible_keys, 1)
                cursors[keys[int(rand_i)]].append(read)

            # Have we passed at least one cursor?
            if read.reference_start > closest_cursor:

                # Check all the cursors
                keys = list(cursors.keys())
                for c_i in np.argsort(keys):
                    c_cursor = keys[c_i]
                    if read.reference_start >= c_cursor:
                        print("%d : passed cursor %d" % (read.reference_start, c_cursor))
                        tracks_at_cursor = np.argwhere(track_ends == c_cursor).flatten()
                        chosen_reads = np.random.choice(cursors[c_cursor], min(len(cursors[c_cursor]), len(tracks_at_cursor)), replace=False)

                        read_i = -1
                        for read_i, chosen_read in enumerate(chosen_reads):
                            n_reads += 1
                            curr_track = tracks_at_cursor[read_i]
                            track_reads[curr_track] = chosen_read
                            track_ends[curr_track] = chosen_read.reference_end 
                            cursors[chosen_read.reference_end] = []
                            track_cov[chosen_read.reference_start - region.start : chosen_read.reference_end - region.start] += 1

                            seen_reads.add(chosen_read.query_name)
                            bam_ret_d["%s:%d" % (chosen_read.query_name, chosen_read.reference_start)] = 1 # send read name to write to new BAM after all threads are done

                            if args.output_fasta:
                                read_ret_q.put(chosen_read.to_string()) # send read data to be written to fasta/fastq

                        #print("BLIP", track_ends)

                        #print("READ_i")
                        #for read_i, read in enumerate(chosen_reads):
                        #    curr_track = tracks_at_cursor[read_i]
                        #    print(read_i, curr_track, track_ends[curr_track])
                        #print("READ_j")
                        #for read_j in range(read_i+1, len(tracks_at_cursor)):
                        #    curr_track = tracks_at_cursor[read_j]
                        #    print(read_j, curr_track, track_ends[curr_track])

                        del cursors[c_cursor] # drop the reads from the cursor watch


                # Count the number of tracks that need to be filled at this position
                closest_cursor = np.amin(track_ends)
                logger.info("\n\nClosest cursor advanced to %d, last read %d" % (closest_cursor, read.reference_start))
                if read.reference_start > closest_cursor:
                    next_i = 1
                    next_cursor = closest_cursor
                    while next_cursor <= read.reference_start:
                        try:
                            next_cursor = sorted(set(track_ends))[next_i]
                            next_i += 1
                        except IndexError:
                            next_cursor = np.inf
                            break

                    next_cursor = min(next_cursor, read.reference_start + 1000)

                    #for read_j in range(read_i+1, len(tracks_at_cursor)):
                    #    curr_track = tracks_at_cursor[read_j]
                    #    track_ends[curr_track] = next_cursor
                    #print("BLOP", track_ends)

                    # prevent assigning a track end that the cursor has already passed
                    for t_i, track_end in enumerate(track_ends):
                        if track_end <= read.reference_start:
                            track_ends[t_i] = next_cursor

                            if next_cursor not in cursors:
                                cursors[next_cursor] = []
                            #cursors[next_cursor].extend(cursors[track_end])
                            #del cursors[track_end]



        median_depth = np.median(track_cov)
        stdv_depth = np.std(track_cov)
        logger.info(u'region: {}, reads: {}, target: {}, depth: {:.0f}X (\u00B1{:.1f}).'.format(region.ref_name, n_reads, args.depth, median_depth, stdv_depth))


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
