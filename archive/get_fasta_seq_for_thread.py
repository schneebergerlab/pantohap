#!/usr/bin/env python3
"""
Given the thread, find the representive genome and fetch its fasta sequence
"""
import argparse
from collections import defaultdict, deque, Counter
from itertools import product
import socket
import string
import pandas as pd
from hometools.hometools import readfasta
from tqdm import tqdm
import os
import sys

class hapobject:
    """
    A haplotype block object
    """

    def __init__(self, id, start, end, genomes):
        self.id = id
        self.start = start
        self.end = end
        self.genomes = genomes

    # END

    def hasgen(self, gen):
        return gen in self.genomes
    # END
# END

def readthreads(f):
    listtoint = lambda x: list(map(int, x.replace('[', '').replace(']', '').split(',')))
    threads = deque()
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split('\t')
            if line[4] == 'NA': continue
            # ts = listtoint(line[2])
            ts = listtoint(line[3].replace('imputed_', ''))
            # Only select threads that cover at least three nodes
            if len(ts) >= 2:
                # threads.append([int(line[1]), listtoint(line[2])])
                threads.append([int(line[2]), ts, line[4].split(), int(line[0])])
    return threads


def get_hap_map_len(df, n):
    """
    df: containing hap lens
    n: node position
    """
    df = df.loc[df.window==n, ['sample_start', 'sample_end']]
    df['len'] = df['sample_end'] - df['sample_start']
    if df['len'].sum() == 0:
        return 0
    else:
        return df['len'].sum() + df.shape[0]
# END

def get_hap_seq(df, windows, fasta, gen, threadid):
    df = df.loc[df.window.isin(windows), ['sample_start', 'sample_end']]
    df = df.loc[~(df.sample_start == 0)]
    s = df.sample_start.min()
    e = df.sample_end.max()
    fseq = readfasta(fasta)
    c = list(fseq.keys())[0]
    seq = list(fseq.values())[0]
    outc = f'{gen}_{c}_{s}_{e}_threadid_{threadid}'
    return outc, seq[(s-1):e]
# END

def main(args):
    parser.add_argument('hap', help='Haplotype graph', type=argparse.FileType('r'))
    parser.add_argument('threads', help='Threads from EM', type=argparse.FileType('r'))
    parser.add_argument('haplen', help='Selected fasta regions for haps', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output output file name', type=argparse.FileType('w'))
    parser.add_argument('chrid', help='Chromosome of interes', type=str)

    HAPFIN = args.hap.name
    THFIN = args.threads.name
    HLFIN = args.haplen.name
    OUTFOUT = args.out.name
    CHRID = args.chrid

    if socket.gethostname() == 'LMBIDP1-LXSCH01':
        # Indir for local compute
        indir = '/home/ra98jam/d16/projects/potato_hap_example/data/'
        pwd = '/home/ra98jam/d16/projects/potato_hap_example/results/threading/thread_fastas/'
    else:
        # Indir for cluster compute
        indir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/'
        pwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/threading/thread_fastas/'
    os.chdir(pwd)
    #
    # HAPFIN = '../../../data/chr10/haplotype_graph_chr10_div0.1_2024_08_02.txt'
    # THFIN = '../threads_forManish_26_08_24/WhiteRose_chr10_div0.1.threads.tsv'
    # HLFIN = '../../kmer_analysis/node_kmers/fasta_len_in_nodes.csv'
    # OUTFOUT = 'tmp.threads.fa'
    # CHRID = 'chr10'

    # hapfin = "/home/ra98jam/d16/projects/potato_hap_example/data/chr06/haplotype_graph_chr06_div0.1_2024_08_02.txt"
    # Read haplotype graph
    hapoblist = deque()
    with open(HAPFIN, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            hapoblist.append(hapobject(int(line[3]), int(line[0]), int(line[1]), sorted(line[2].split(','))))

    # thfin = '/home/ra98jam/d16/projects/potato_hap_example/results/threading/threads_forManish_26_08_24/RussetBurbank_chr06_div0.1.threads.tsv'
    threads = readthreads(THFIN)

    # hlfin = "/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/node_kmers/fasta_len_in_nodes.csv" # File containing the coordinates of mapped regions
    hliter = pd.read_table(HLFIN, iterator=True, chunksize=10000)
    hldf = pd.concat([chunk[chunk['chromosome'] == CHRID] for chunk in hliter]).groupby('sample')

    # samples = [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(1, 5))]

    maxlen = deque()
    for thread in threads:
        slens = defaultdict(int)
        tnodes = thread[1]
        nodegen = {t: hapoblist[t].genomes for t in tnodes}
        for k, v in nodegen.items():
            for gen in v:
                if gen == 'dm': continue
                # print(k, gen, get_hap_map_len(hldf.get_group(gen), hapoblist[k].start))
                slens[gen] += get_hap_map_len(hldf.get_group(gen), hapoblist[k].start)

        candidates = thread[2][0].split(',')
        badcan = [t for t in candidates if slens[t] == 0]
        if len(badcan) == len(candidates):
            maxlen.append(sorted(slens, key=lambda x: slens[x])[-1])
            continue

        while True:
            maxlenhap = sorted(slens, key=lambda x: slens[x])[-1]
            if maxlenhap in thread[2][0]:
                # print(maxlenhap)
                maxlen.append(maxlenhap)
                break
            else:
                slens.pop(maxlenhap)
            if len(slens) == 0:
                break
    maxlen = list(maxlen)

    # for i, thread in enumerate(threads):
    #     try:
    #         assert maxlen[i] in thread[2][0]
    #     except AssertionError:
    #         print(i, threads[i], maxlen[i])


    with open(OUTFOUT, 'w') as fout:
        for i, thread in tqdm(enumerate(threads)):
            threadid = thread[3]
            gen = maxlen[i]
            s, h = gen.split('_hap')
            windows = [hapoblist[n].start for n in thread[1]]
            try:
                outc, outseq = get_hap_seq(hldf.get_group(gen), windows, f'{indir}/{CHRID}/{s}_{CHRID}_hap{h}.fa', gen, threadid)
            except:
                print(thread, HAPFIN, THFIN)
                sys.exit()
            outseq = "\n".join([outseq[i:i+60] for i in range(0, len(outseq), 60)])
            fout.write('>{outc}\n{outseq}\n'.format(outc=outc, outseq=outseq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get pseudo assembly contigs')
    parser.add_argument('hap', help='Haplotype graph', type=argparse.FileType('r'))
    parser.add_argument('threads', help='Threads from EM', type=argparse.FileType('r'))
    parser.add_argument('haplen', help='Selected fasta regions for haps', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output output file name', type=argparse.FileType('w'))
    parser.add_argument('chrid', help='Chromosome of interes', type=str)

    args = parser.parse_args()
    main(args)