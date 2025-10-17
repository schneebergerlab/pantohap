#!/usr/bin/env python3

"""
For a given sample, fetch kmer counts for node-specific and window-specific kmers.
"""
import argparse

def isgzip(f):
    '''
    Checks if a given file is gzipped. Returns true if the is gzipped.
    '''
    from gzip import open as gzopen
    from gzip import BadGzipFile
    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgz = True
        except BadGzipFile:
            isgz = False
    return isgz
# END


def count_kmers_from_samples(samplekmersfin):
    """
    For a given sample, fetch kmer counts for node-specific and window-specific kmers.
    :return: Generates
    """
    import os
    # from hometools.hometools import isgzip
    from gzip import open as gzopen
    from collections import defaultdict, deque
    import numpy as np
    import pandas as pd
    from tqdm import tqdm
    import logging

    logger = logging.getLogger(__name__)
    samplekmersfin = args.sfin.name
    nodekmersfin = args.kfin.name
    outfin = args.out.name if args.out is not None else 'node_k_stats.txt'
    logger.info(f'Writing output to {outfin}')

    # Check if gzip
    fo = gzopen if isgzip(samplekmersfin) else open

    with open(nodekmersfin, 'r') as fin:
        kmers = {line.strip().split()[1].encode(): 0 for line in fin}

    # Get kmer counts for sample
    with fo(samplekmersfin, 'r') as fin:
        for i, line in tqdm(enumerate(fin)):
            line = line.strip().split()
            if line[0] in kmers:
                kmers[line[0]] = int(line[1])

    absentkmer = [k for k in kmers if kmers[k] == 0]
    for kmer in absentkmer:
        kmers.pop(kmer)
    del absentkmer

    # for each node get distribution of kmer counts
    nodeks = defaultdict(int)           # Number of marker kmers in node
    nodekmercnts = defaultdict(deque)   # Number of reads with the specific kmer
    with open(nodekmersfin, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            nodeks[(int(line[0]))] += 1
            try:
                nodekmercnts[int(line[0])].append(kmers[line[1].encode()])
            except KeyError:
                pass
            except Exception as e:
                if line == []:
                    pass
                else:
                    raise e
    with open(outfin, 'w') as fout:
        for k in nodekmercnts.keys():
            m = round(np.mean(nodekmercnts[k]), 2) if len(nodekmercnts[k]) > 0 else 0
            fout.write(f'{k}\t{nodeks[k]}\t{len(nodekmercnts[k])}\t{round(len(nodekmercnts[k])/nodeks[k], 2)}\t{m}\t{",".join(list(map(str, nodekmercnts[k])))}\n')

    nodestats = pd.read_table(outfin, header=None)
    nodestats.sort_values([2, 3], ascending=False, inplace=True)
    nodestats.to_csv(outfin, sep='\t', header=False, index=False)

    logger.info('Finished')
    # Generate plots here: analysis_plots.node_k_stats_plots
    return
# END


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get Kmer stats for samples.')
    parser.add_argument(help='Sample kmer counts', dest='sfin', type=argparse.FileType('r'))
    parser.add_argument(help='List of node kmers', dest='kfin', type=argparse.FileType('r'))
    parser.add_argument('--out', dest='out', help='Output file name', type=argparse.FileType('w'))
    args = parser.parse_args()
    count_kmers_from_samples(args)