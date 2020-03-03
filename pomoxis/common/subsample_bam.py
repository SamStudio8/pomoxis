import argparse
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Process, Queue, Array

import functools
import logging
import multiprocessing
import os

from intervaltree import IntervalTree, Interval
import numpy as np
import pysam
from random import sample

from pomoxis.common.util import parse_regions, Region
from pomoxis.common.coverage_from_bam import coverage_summary_of_region
from pomoxis.common.stats_from_bam import stats_from_aligned_read


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('subsample bam to uniform or proportional depth',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam',
        help='input bam file.')
    parser.add_argument('depth', nargs='+', type=int,
        help='Target depth.')
    parser.add_argument('-o', '--output_prefix', default='sub_sampled',
        help='Output prefix')
    parser.add_argument('--output-override', default=None,
        help='Specify a path to override generation of one BAM per region')
    parser.add_argument('--output-fastq', default=None,
        help='Create an additional file containing the query_name of subsampled reads.\nRequires --output-override.')
    parser.add_argument('--output-fasta', default=None,
        help='Create an additional file containing the query_name of subsampled reads.\nRequires --output-override.')

    iparser = parser.add_mutually_exclusive_group()
    jparser = iparser.add_mutually_exclusive_group()
    jparser.add_argument('-t', '--threads', type=int, default=-1,
        help='Number of threads to use.')
    #jparser.add_argument('--no-index', action='store_true',
    #    help='Stream SAM or BAM input without index (e.g. via stdin), only supports one thread')
    iparser.add_argument('-r', '--regions', nargs='+',
        help='Only process given regions (requires index).')

    parser.add_argument('-p', '--profile', type=int, default=1000,
        help='Stride in genomic coordinates for depth profile.')
    parser.add_argument('-O', '--orientation', choices=['fwd', 'rev'],
        help='Sample only forward or reverse reads.')
    parser.add_argument('-q', '--quality', type=float,
        help='Filter reads by mean qscore.')
    parser.add_argument('-a', '--accuracy', type=float,
        help='Filter reads by accuracy.')
    parser.add_argument('-c', '--coverage', type=float,
        help='Filter reads by coverage (what fraction of the read aligns).')

    eparser = parser.add_mutually_exclusive_group()
    eparser.add_argument('--any_fail', action='store_true',
        help='Exit with an error if any region has insufficient coverage.')
    eparser.add_argument('--all_fail', action='store_true',
        help='Exit with an error if all regions have insufficient coverage.')
    eparser.add_argument('--output-anyway', action='store_true',
        help='Output a region even if the reads were of insufficient depth.')

    uparser = parser.add_argument_group('Uniform sampling options')
    uparser.add_argument('-x', '--patience', default=5, type=int,
        help='Maximum iterations with no change in median coverage before aborting.')
    uparser.add_argument('-s', '--stride', type=int, default=1000,
        help='Stride in genomic coordinates when searching for new reads. Smaller can lead to more compact pileup.')

    pparser = parser.add_argument_group('Proportional sampling options')
    pparser.add_argument('-P', '--proportional', default=False, action='store_true',
        help='Activate proportional sampling, thus keeping depth variations of the pileup.')
    pparser.add_argument('-S', '--seed', default=None, type=int,
        help='Random seed for proportional downsampling of reads.')

    args = parser.parse_args()
    if args.threads == -1:
        args.threads = multiprocessing.cpu_count()

    with pysam.AlignmentFile(args.bam) as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))

        if args.regions is not None:
            regions = parse_regions(args.regions, ref_lengths=ref_lengths)
        else:
            regions = [Region(ref_name=r, start=0, end=ref_lengths[r]) for r in bam.references]

    #if args.proportional:
    #    worker = functools.partial(subsample_region_proportionally, args=args)
    #else:
    #    worker = functools.partial(subsample_region_uniformly, args=args)

    src_bam = pysam.AlignmentFile(args.bam, "rb")
    out_bam = pysam.AlignmentFile('-', "wh", template=src_bam)
    src_bam.close()

    enough_depth = []
    reads = []
    #with ProcessPoolExecutor(max_workers=args.threads) as executor:
    #    for res in executor.map(worker, regions):
    #        enough_depth.append(res[0])
    #        reads.extend(res[1])
    work_queue = Queue() 
    processes = []

    #manager = multiprocessing.Manager()
    #return_queue = manager.dict()
    return_queue = Queue()
    for _ in range(args.threads):
        p = Process(target=subsample_region_uniformly, args=(work_queue,return_queue,args))
        processes.append(p)
    p = Process(target=spew_reads, args=(return_queue,args))
    processes.append(p)

    for p in processes:
        p.start()

    for region in regions:
        work_queue.put({"region": region})

    for _ in range(args.threads):
        work_queue.put(None)

    # Wait for processes to complete work
    for p in processes:
        p.join()

    logger = logging.getLogger()
    logger.info("Pulling reads back together")
    out_bam.close()

    #if args.output_override:
    #    _write_bam(args.bam, '', regions, return_queue, override=args.output_override, override_fastq=args.output_fastq)
    #    #_write_coverage(prefix, '', '', args.profile)

    if args.any_fail and not all(enough_depth):
        raise RuntimeError('Insufficient read coverage for one or more requested regions.')
    if args.all_fail and all([not x for x in enough_depth]):
        raise RuntimeError('Insufficient read coverage for all requested regions.')

def spew_reads(return_q, args):
    if args.output_fastq:
        ofq_fh = open(args.output_fastq, 'w')
    if args.output_fasta:
        ofa_fh = open(args.output_fasta, 'w')
    written = 0
    dead = 0

    while True:
        read = return_q.get()
        if read is None:
            dead += 1
            if dead == args.threads:
                return None
            continue

        written += 1
        if args.output_fastq:
            fields = read.split()
            ofq_fh.write("@%s\n%s\n+\n%s\n" % (fields[0], read.split()[9], read.split()[10]))
        if args.output_fasta:
            fields = read.split()
            ofa_fh.write(">%s\n%s\n" % (fields[0], read.split()[9]))
        print(read)

    if args.output_fastq:
        ofq_fh.close()
    if args.output_fasta:
        ofa_fh.close()


def filter_read(r, args):
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
            return True

    # don't filter
    return False


def subsample_region_uniformly(work_q, return_q, args):
    while True:
        work = work_q.get()
        if work is None:
            return_q.put(None)
            return
        else:
            region = work["region"]

        logger = logging.getLogger(region.ref_name)
        logger.info("PULLED %s FROM QUEUE" % region.ref_name)

        bam = pysam.AlignmentFile(args.bam)
        coverage = np.zeros(region.end - region.start, dtype=np.uint16)

        # begin
        #tracks = []
        reads = 0
        track_ends = {}
        #pileups = bam.pileup(contig=region.ref_name)
        #pileupcolumn = next(pileups); cursor = 0
        #TODO I want to bring this back to get a better starting read set

        #open_reads = [x for x in pileupcolumn.pileups if (x.alignment.reference_start < 100)]
        #for i,r in enumerate(sample(open_reads, args.depth[0])):
        #    tracks.append( [r.alignment] )
        #    track_ends[i] = r.alignment.reference_end
        #    coverage[r.alignment.reference_start : r.alignment.reference_end] += 1
        for i in range(args.depth[0]):
            #tracks.append( [] )
            track_ends[i] = 0
        cursor = 0

        #return_q[region.ref_name] = {}
        for read in bam.fetch(contig=region.ref_name):
            # TODO We want to get some randomness in here
            if read.reference_start > cursor:
                if filter_read(read, args):
                    continue
                this_track = [i for i in track_ends if track_ends[i] == cursor][0] # select one if they have same start
                #tracks[this_track].append(read.query_name+'@'+str(read.reference_start))
                #return_q[region.ref_name][read.query_name+'@'+str(read.reference_start)] = 1
                return_q.put(read.to_string())
                #reads.append(str(read))
                track_ends[this_track] = read.reference_end
                #out_bam.write(read)
                coverage[read.reference_start : read.reference_end] += 1
                reads += 1
                #for t in this_tracks:
                #    for i, r in enumerate(sample(tree.at(cursor), len(this_tracks))):
                #        tracks[this_tracks[i]].append(r)
                #        track_ends[this_tracks[i]] = r.reference_end
                #        #out_bam.write(r.alignment)
                #        coverage[r.reference_start : r.reference_end] += 1
                cursor = min(track_ends.values()) 

        #reads = []
        #[reads.extend(track) for track in tracks]
        #return_q.extend(reads)

        median_depth = np.median(coverage)
        stdv_depth = np.std(coverage)
        logger.info(u'reads: {}, depth: {:.0f}X (\u00B1{:.1f}).'.format(
            reads, median_depth, stdv_depth))

        logger.info("PUSHED %s TO QUEUE" % region.ref_name)


def _nearest_overlapping_point(src, point, step, lengths={}):
    """Find the interval with the closest start point to a given point.

    :param src: IntervalTree instance.
    :param point: query point.

    :returns: Interval instance of interval with closest start.

    """
    items = src.at(point)
    #items = src.envelop(point, point+step)
    if len(items) > 0:
        return next(iter(items))

    #if len(items) == 0:
    #    return None
    #
    #if lengths:
    #    pass
    #    #items = sorted(items, key=lambda x: lengths[x.data], reverse=True) # this actually codes a preference for longest reads first, which we might not want
    #else:
    #    items = sorted(items, key=lambda x: x.end - x.begin, reverse=True)
    #items.sort(key=lambda x: abs(x.begin - point))
    #return items[0]


def _write_bam(bam, prefix, region, sequences, override=False, override_fastq=None):
    # filtered bam

    logger = logging.getLogger()
    if override:
        if override_fastq:
            ofq_fh = open(override_fastq, 'w')

        logger.info("Writing BAM.")
        written = 0
        output = override
        src_bam = pysam.AlignmentFile(bam, "rb")
        out_bam = pysam.AlignmentFile(output, "wb", template=src_bam)
        for r in region:
            for read in src_bam.fetch(contig=r.ref_name):
                if '%s@%d' % (read.query_name, read.reference_start) in sequences[read.reference_name]:
                    out_bam.write(read)
                    written+=1
                    if override_fastq:
                        ofq_fh.write("%s\n" % read.query_name)

                    if written % 10000 == 0:
                        logger.info("%.2f written..." % (written/len(sequences)))

        if override_fastq:
            ofq_fh.close()
    else:
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
    if override != '-':
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
