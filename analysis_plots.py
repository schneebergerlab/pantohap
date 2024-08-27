import matplotlib.pyplot as plt


def summary_plots_kmers_per_node():
    """
    Collection of plots to describe the distribution of kmers in nodes, windows, and chromosomes
    :return:
    """
    from collections import deque, Counter
    from matplotlib import pyplot as plt
    from tqdm import tqdm
    from hometools.plot import cleanax

    allnodes = tuple(range(2948))    # Number of nodes as identified in /home/ra98jam/d16/projects/potato_hap_example/data/haplotype_graph.txt
    # <editor-fold desc="Kmer count distribution for the nodes">
    for k in (21, 31, 41, 51):
        nodecnts = deque()
        with open('/home/ra98jam/projects/pantohap/results/potato_haps/nodekmers.txt', 'r') as fin:
            for line in tqdm(fin):
                try:
                    nodecnts.append(int(line.strip().split()[0]))
                except IndexError:
                    pass
        nodecnts = Counter(nodecnts)
        nodecntscomp = defaultdict(int)
        for n in allnodes:
            try:
                nodecntscomp[n] = nodecnts[n]
            except KeyError:
                pass
        # kmer Count histogram
        plt.hist(nodecntscomp.values(), bins=200)
        plt.yscale('log')
        plt.xscale('linear')
        plt.xlabel('Number of Kmers (k=51)')
        plt.ylabel('Number of Nodes')
        plt.title('Total number of nodes in Chr02')
        plt.tight_layout()
        plt.savefig('/home/ra98jam/projects/pantohap/results/potato_haps/kmer_frequencies_per_node.pdf')
        plt.close()

    # </editor-fold>

    # Kmer distribution along the chromosome
    ## Read node location
    hapg = pd.read_table(f'/home/ra98jam/projects/pantohap/results/potato_haps/haplotype_graph.txt', header=None)
    kmerfreq = [10, 100, 500, 1000, 5000, 10000]
    fig = plt.figure(figsize=[12, 8])
    addleg = True
    for i, k in enumerate(kmerfreq):
        chrkmercnt = dict()
        for grp in hapg.groupby([0, 1]):
            ncnt = grp[1].shape[0]
            highkmernct = sum([True for n in grp[1][3] if nodecntscomp[n] > k])
            chrkmercnt[grp[0]] = (ncnt, highkmernct)
        xranges = list(chrkmercnt.keys())
        ax = fig.add_subplot(len(kmerfreq), 1, i+1)
        ax.plot([np.mean(x) for x in xranges], [chrkmercnt[x][0] for x in xranges], color='black', label='In graph')
        ax.plot([np.mean(x) for x in xranges], [chrkmercnt[x][1] for x in xranges], color='red', label='High Kmer count')
        ax.set_xlabel('Chromosome position')
        ax.set_ylabel('#Nodes')
        ax.set_title(f'Kmer freq cutoff {k}, k=51')
        ax = cleanax(ax)
        if addleg:
            ax.legend()
            addleg = False
    plt.tight_layout()
    plt.savefig('/home/ra98jam/projects/pantohap/results/potato_haps/kmer_cutoff_along_chromosome.pdf')
    plt.close()
    return
# END

def node_fasta_sequence_stats():
    from collections import defaultdict, deque
    import socket
    import os
    from datetime import datetime
    import seaborn as sns
    from multiprocessing import Pool
    from itertools import product
    import string
    import numpy as np
    from hometools.hometools import unlist
    import pandas as pd

    # Set paths
    if socket.gethostname() == 'LMBIDP1-LXSCH01':
        # Indir for local compute
        indir = '/home/ra98jam/d16/projects/potato_hap_example/data/'
        pwd = '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/node_kmers/'
    else:
        # Indir for cluster compute
        indir = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/'
        pwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/node_kmers/'
    # set variables
    chrs = ["chr{:02d}".format(i) for i in range(1, 13)]
    samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(1, 5))]  # Without Otava
    WINDOW_SIZE = 100000

    # Read chromosome lengths
    chrlengths = deque()
    for c in chrs:
        for s in string.ascii_uppercase[:10]:
            for h in range(1, 5):
                _, l, _, _, _ = np.loadtxt(f'{indir}/{c}/{s}_{c}_hap{h}.fa.fai', dtype=str)
                chrlengths.append((c, f'{s}_hap{h}', int(l)))
    dmlengths = defaultdict()
    for c in chrs:
        _, l, _, _, _ = np.loadtxt(f'{indir}/{c}/DM_{c}.fa.fai', dtype=str)
        dmlengths[c] = int(l)

    # Get length of syntenic regions
    def get_fasta_len(c):
        fasta_len = deque()
        for s in samples[1:3]:
            print(c, s, str(datetime.now()))
            windows = np.arange(1, dmlengths[c], WINDOW_SIZE)
            indir = f'{pwd}/{c}/{s}/synfastas/'
            for w in windows:
                f = f"syn_fasta_{s}_bin_{w}.fasta"
                try:
                    with open(f'{indir}/{f}', 'r') as fin:
                        for line in fin:
                            if line[0] == '>':
                                line = line.strip().split('_')
                                fasta_len.append([c, s, w] + [line[4], line[5], line[7], line[8]])
                except FileNotFoundError:
                    fasta_len.append([c, s, w] + [0, 0, 0, 0])
        return fasta_len

    with Pool(processes=12) as pool:
        fasta_len2 = pool.map(get_fasta_len, chrs)

    tmp = [j for i in fasta_len2 for j in i]
    fasta_len = tmp

    fasta_len_df = pd.DataFrame(fasta_len)
    fasta_len_df.columns = 'chromosome sample window dm_start dm_end sample_start sample_end'.split()
    fasta_len_df['window dm_start dm_end sample_start sample_end'.split()] = fasta_len_df['window dm_start dm_end sample_start sample_end'.split()].astype(int)
    fasta_len_df.sort_values('chromosome sample window dm_start'.split(), ascending=True, inplace=True)
    fasta_len_df.to_csv(f'{pwd}/fasta_len_in_nodes.csv', index=False, sep='\t')

    # Get fasta len summary statistics
    fasta_len_df['dm_len'] = (fasta_len_df['dm_end'] - fasta_len_df['dm_start'] + 1)/WINDOW_SIZE
    fasta_len_df['sample_len'] = (fasta_len_df['sample_end'] - fasta_len_df['sample_start'] + 1)/WINDOW_SIZE
    flagg = fasta_len_df.groupby('chromosome sample window'.split()).agg({'dm_len': 'sum', 'sample_len': 'sum'})
    for c in chrs:
        pltdf = flagg.loc[(c)]
        pltdf.reset_index(inplace=True)
        fig, ax = plt.subplots(figsize=(10, 3))
        ax = sns.lineplot(pltdf, x='window', y='sample_len', ax=ax)
        ax = sns.scatterplot(pltdf, x='window', y='sample_len', s=2, ax=ax)
        ax.set_title(c)
        ax.set_ylim([-0.1, 1.1])
        plt.tight_layout()
        plt.savefig(f'{pwd}/{c}_fasta_len_dist.png')
        plt.close()

    fasta_len_df['dm_len'] = fasta_len_df['dm_end'] - fasta_len_df['dm_start'] + 1
    fasta_len_df['sample_len'] = fasta_len_df['sample_end'] - fasta_len_df['sample_start'] + 1
    flagg2 = fasta_len_df.groupby('chromosome sample'.split()).agg({'dm_len': 'sum', 'sample_len': 'sum'})
    # flagg2.reset_index(inplace=True)
    chrldf = pd.DataFrame(chrlengths)
    chrldf.columns = ['chromosome', 'sample', 'chr_len']
    chrldf.set_index(['chromosome', 'sample'], inplace=True)
    flagg2.sample_len = flagg2.sample_len / chrldf.chr_len
    flagg2.dm_len = flagg2.dm_len / chrldf.chr_len
    flagg2.reset_index(level=1, inplace=True)
    flagg2['Genome'] = flagg2['sample'].str.split('_').str[0]
    flagg2['Haplotype'] = flagg2['sample'].str.split('_').str[1]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax = sns.stripplot(flagg2, x='chromosome', y='sample_len', jitter=0.25, hue='Genome', ax=ax)
    ax.set_ylim([-0.1, 1.1])
    plt.tight_layout()
    plt.savefig(f'{pwd}/all_chr_region_covered.png')
    plt.close()
    return
# END

def node_k_stats_plots():

    nodestats = pd.read_table('/home/ra98jam/projects/pantohap/results/potato_haps/node_k_stats.txt', header=None)
    nodestats.sort_values([2, 3], ascending=False, inplace=True)
    nodestats.to_csv('/home/ra98jam/projects/pantohap/results/potato_haps/node_k_stats.txt', sep='\t', header=False, index=False)

    fig = sns.jointplot(nodestats, x=1, y=3)
    plt.xlabel('Number of 51-mers supporting node')
    plt.ylabel('Percent of 51-mers identified')
    plt.tight_layout()
    plt.savefig('/home/ra98jam/projects/pantohap/results/potato_haps/total_kmer_counts_vs_ratio.pdf')

    pltdf = nodestats.copy()
    pltdf.columns = 'k kmer_total kmer_found kmer_ratio avg_kmer_cnt kmer_cnts'.split()
    pltdf = pltdf.loc[pltdf.kmer_ratio > 0.8]
    nodepos = pd.read_table('/home/ra98jam/projects/pantohap/results/potato_haps/haplotype_graph.txt', header=None)
    nodepos.columns = 'start end haps k'.split()
    nodepos['xpos'] = (nodepos['start'] + nodepos['end'])/2
    pltdf = pltdf.merge(nodepos, on='k', how='inner')
    ypos = deque()
    for grp in pltdf.groupby('xpos'):
        ypos.extend(list(range(grp[1].shape[0])))
    pltdf['ypos'] = ypos

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)
    ax1 = sns.scatterplot(pltdf, x='xpos', y='ypos', hue='kmer_found', ax=ax1)
    ax2 = sns.scatterplot(pltdf, x='xpos', y='ypos', hue='kmer_ratio', ax=ax2)
    ax3 = sns.scatterplot(pltdf, x='xpos', y='ypos', hue='avg_kmer_cnt', ax=ax3)
    ax1.set_title('Number of Kmers found in Otava')
    ax2.set_title('Proportion of unique found in Otava')
    ax3.set_title('Average K-mer count in Otava reads')
    plt.tight_layout()
    plt.savefig('/home/ra98jam/projects/pantohap/results/potato_haps/kmer_dist_along_chromosome.pdf')

    return
# END


def plot_haplotype_graph():
    """
    Plot a given a haplotype graph
    :return:
    """
    import igraph as ig
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

    hapfin = '/home/ra98jam/projects/pantohap/results/potato_haps/haplotype_graph.txt'
    hapoblist = deque()
    with open(hapfin, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            hapoblist.append(hapobject(int(line[3]), int(line[0]), int(line[1]), sorted(line[2].split(','))))
    #
    # cnt = 0
    # for i, h in enumerate(hapblocks):
    #     for ss in h[2:]:
    #         ss = sorted(ss.split(','))
    #         hapoblist.append(hapobject(cnt, h[0], h[1] - 1, ss))
    #         cnt += 1

    # Create a Graph object for the haplotype blocks and use it to determine the Y-coordinate for visualisation of the blocks
    G = ig.Graph()
    G = ig.Graph.as_directed(G)
    G.add_vertices(range(len(hapoblist)))
    addededge = set()
    startnodes = deque()
    for hap in samples:
        hc = -1
        for h in hapoblist:
            if h.hasgen(hap):
                if hc == -1:
                    hc = h.id
                    startnodes.append(h.id)
                    continue
                if (hc, h.id) not in addededge:
                    G.add_edge(hc, h.id)
                    addededge.add((hc, h.id))
                hc = h.id
    # to_delete_ids = [v.index for v in G.vs if v.degree() == 0]
    # G.delete_vertices(to_delete_ids)
    G.vs['dist'] = -1
    G.vs['height'] = -1
    G.vs['start'] = [h.start for h in hapoblist]

    # roots = [v.index for v in G.vs if v['dist'] == 0]
    roots = [v.index for v in G.vs if v['start'] == 1]
    G.vs[roots]['height'] = np.ceil(np.arange(len(roots)) - (len(roots) / 2))
    # maxx = max(G.vs['dist'])
    # for i in range(1, maxx):
    for i in sorted(set(G.vs['start']))[1:]:
        ns = [v.index for v in G.vs if v['start'] == i]
        heights = np.ceil(np.arange(len(ns)) - (len(ns) / 2))
        hord = sorted(range(len(ns)), key=lambda x: np.mean(G.vs[G.predecessors(ns[x])]['height']))
        G.vs[ns]['height'] = [heights[hord.index(j)] for j in range(len(ns))]

    G.vs['height'] += np.random.uniform(-0.25, 0.25, len(G.vs))

    fig, ax = plt.subplots(figsize=[8, 6])
    plt.tight_layout()
    ax = cleanax(ax)
    ax.set_xlim(0, chrsize)
    # ax.set_xlim(34800000, 35400000)

    segs = deque()
    cs = deque()
    for i, h in enumerate(list(hapoblist)):
        segs.append((((h.start + 10000), G.vs[i]['height']), ((h.end - 10000), G.vs[i]['height'])))
        if 'dm' in h.genomes:
            cs.append('lightgrey')
        else:
            cs.append('black')

    line_segments = LineCollection(segs, colors=cs, linestyle='solid', linewidths=1)
    ax.add_collection(line_segments)
    ax.set_ylim(-5, 7)


    # plt.savefig('potato_haplotypes.png')
    # starts = list(range(34500000, 37500000, 100000))
    starts = list(range(34800000, 35400000, 100000))

    return
# END

