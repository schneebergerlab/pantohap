"""
Given the thread, find the representive genome and fetch its fasta sequence
"""
from collections import defaultdict, deque, Counter
from itertools import product
import socket
import string
import pandas as pd
from hometools.hometools import readfasta
from tqdm import tqdm
 igraph as ig
import math
import numpy as np
from numpy.polynomial import Polynomial as P
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib import colors as mcol
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
from matplotlib.pyplot import get_cmap
from hometools.hometools import mylogger
from hometools.plot import cleanax
import os
from hometools.plot import plotdensity
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
    df['len'] = df['sample_end'] - df['sample_start'] + 1
    return sum(df['len'])
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


hapfin = "/home/ra98jam/d16/projects/potato_hap_example/data/chr06/haplotype_graph_chr06_div0.1_2024_08_02.txt"
# Read haplotype graph
hapoblist = deque()
with open(hapfin, 'r') as fin:
    for line in fin:
        line = line.strip().split()
        hapoblist.append(hapobject(int(line[3]), int(line[0]), int(line[1]), sorted(line[2].split(','))))

thfin = '/home/ra98jam/d16/projects/potato_hap_example/results/threading/threads_forManish_26_08_24/RussetBurbank_chr06_div0.1.threads.tsv'
threads = readthreads(thfin)

hlfin = "/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/node_kmers/fasta_len_in_nodes.csv" # File containing the coordinates of mapped regions
hliter = pd.read_table(hlfin, iterator=True, chunksize=10000)
hldf = pd.concat([chunk[chunk['chromosome'] == 'chr06'] for chunk in hliter]).groupby('sample')

samples = [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(1, 5))]

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

for i, thread in enumerate(threads):
    try:
        assert maxlen[i] in thread[2][0]
    except AssertionError:
        print(i, threads[i], maxlen[i])

if socket.gethostname() == 'LMBIDP1-LXSCH01':
    # Indir for local compute
    indir = '/home/ra98jam/d16/projects/potato_hap_example/data/'
    pwd = '/home/ra98jam/d16/projects/potato_hap_example/results/threading/'
else:
    # Indir for cluster compute
    indir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/'
    pwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/threading/'

os.chdir(pwd)
curchr = 'chr06'
with open('tmp.threads.fa', 'w') as fout:
    for i, thread in tqdm(enumerate(threads)):
        threadid = thread[3]
        gen = maxlen[i]
        s, h = gen.split('_hap')
        windows = [hapoblist[n].start for n in thread[1]]
        outc, outseq = get_hap_seq(hldf.get_group(gen), windows, f'{indir}/{curchr}/{s}_{curchr}_hap{h}.fa', gen, threadid)
        outseq = "\n".join([outseq[i:i+60] for i in range(0, len(outseq), 60)])
        fout.write('>{outc}\n{outseq}\n'.format(outc=outc, outseq=outseq))


