# <editor-fold desc='Define imports'>
import string
from collections import defaultdict, deque
from itertools import product
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import igraph as ig
from gzip import open as gzopen
from tqdm import tqdm
import seaborn as sns

from hometools.plot import cleanax, loghist
from hometools.hometools import isgzip
# </editor-fold>


# <editor-fold desc='Define constants'>
chrsize = 46102915
samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]
# </editor-fold>


# <editor-fold desc='Classes'>
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

# </editor-fold>

def editdist(iter1, iter2):
    """
    Get the number of positions at which the two iterators are different.
    Requires that the iterators are of the same length.
    :param str1:
    :param str2:    :return:
    """
    try:
        assert len(iter1) == len(iter2)
    except AssertionError:
        raise AssertionError('The two iterators are of different length.')
    return sum([i!=j for i, j in zip(iter1, iter2)])
# END

def hapnodesfromvcffile(snpfin):
    import pandas as pd
    import numpy as np
    from collections import deque
    from itertools import product
    from tqdm import tqdm
    import igraph as ig
    from hometools.hometools import extractseq, readfasta_iterator, mapbp
    import string
    import pyranges as pr

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


    snpfin = f"/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt"
    # snpfin = f"dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt"
    snpdf = pd.read_table(snpfin)
    snpdf.columns = 'chr position reference alternate'.split() + list(snpdf.columns[4:])
    snpdf['dm'] = '.'
    snpdf.replace('.', 0, inplace=True)
    snpdf.replace('1', 1, inplace=True)

    chrsize = 46102915
    samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]

    a = np.arange(1, chrsize, 100000)
    a = np.append(a, chrsize)

    hapblocks = deque()
    for i, start in tqdm(enumerate(a[:-1])):
        end = a[i + 1]
        snps = snpdf.loc[(snpdf.position >= start) & (snpdf.position < end)].copy()
        snpgt = snps.iloc[:, 7:]
        snpgt = snpgt.T
        if snpgt.shape[1] == 0:
            continue

        haplist = [grp for i, grp in snpgt.groupby(list(snpgt.columns))]
        nhap = len(haplist)

        # Find groups of haps that are very similar (divergence < 10%)
        dist = [editdist(i.iloc[0], j.iloc[0]) for i, j in product(haplist, haplist)]
        dmat = np.reshape([editdist(i.iloc[0], j.iloc[0]) for i, j in product(haplist, haplist)], [nhap, nhap])
        conn = (dmat / haplist[0].shape[1]) < 0.10
        conn = np.triu(conn, 1)
        g = ig.Graph.Adjacency(conn)
        hapgrp = g.connected_components(mode='weak')

        # Merge haps that are within the same group
        hapglist = []
        for hg in hapgrp:
            hap = pd.concat([haplist[i] for i in hg])
            hapglist.append(hap)
        nhap = len(hapglist)
        hapsample = [','.join(h.index.values) for h in hapglist]
        hapblocks.append([start, end] + hapsample)

    # hapoblist = deque()
    # cnt = 0
    # for i, h in enumerate(hapblocks):
    #     for ss in h[2:]:
    #         ss = sorted(ss.split(','))
    #         hapoblist.append(hapobject(cnt, h[0], h[1] - 1, ss))
    #         cnt += 1
    index = 0
    with open(f"/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/haplotype_graph.txt", 'w') as fout:
        for hapblock in hapblocks:
            for hap in hapblock[2:]:
                fout.write(f'{hapblock[0]}\t{hapblock[1]}\t{hap}\t{index}\n')
                index += 1

    return
# END


def plot_haplotype_graph():
    """
    Plot a given a haplotype graph
    :return:
    """
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


def get_node_query_sequence():
    """
    Get query sequence fasta file corresponding to 100kb windows in the reference genome
    :return:
    """
    import pandas as pd
    import numpy as np
    from collections import deque
    from itertools import product
    from tqdm import tqdm
    import igraph as ig
    from hometools.hometools import extractseq, readfasta_iterator, mapbp
    import string
    import pyranges as pr

    chrsize = 46102915
    samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]

    a = np.arange(1, chrsize, 100000)
    a = np.append(a, chrsize)

    # Run in /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/
    # For each sample, select syntenic region in each window (100kb). Skip DM.
    for sample in samples[1:]:
        print(sample)
        g, h = sample.split('_hap')
        # get syn regions
        synr = pd.read_csv(f'dm_{g}_chr02_hap{h}syri.out', header=None, sep='\t')
        synr = synr.loc[synr[10] == 'SYNAL']
        synr[[1, 2, 6, 7]] = synr[[1, 2, 6, 7]].astype(int)
        synd = dict()
        for i in a:
            sinr = synr.loc[(synr[1] < (i+100000)) & (synr[2] > i)]
            if sinr.shape[0] > 0:
                pr1 = pr.from_dict({"Chromosome": ["chr02"],
                                    "Start": [i],
                                    "End": [i+100000]})
                pr2 = pr.from_dict({"Chromosome": ["chr02"]*sinr.shape[0],
                                    "Start": sinr[1],
                                    "End": sinr[2],
                                    "query": sinr[[6, 7]].values.tolist()})
                synd[i] = pr2.intersect(pr1)

        (chrom, seq) = next(readfasta_iterator(open(f'{g}_chr02_hap{h}.fa', 'r')))

        for i, r in synd.items():
            with open(f'/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/synfastas/syn_fasta_{sample}_bin_{i}.fasta',
                    'w') as fout:
                for c in r.df.itertuples(index=False):
                    fout.write(f'>chr02_{g}_hap{h}_r_{c[1]}_{c[2]}_q_{c[3][0]}_{c[3][1]}\n')
                    fout.write(seq[(c[3][0]-1):c[3][1]] + '\n')
    return
# END


def get_unique_kmers_per_node(cwd, k):
    """
    Read the haplotype graph and:
    1) Get set of kmers for each node
    2) Filter out kmers that are present in more than 1 node
    3) Save the unique kmers for each node
    :return:
    """

    def gethapmers(row):
        b = row[0]
        haps = row[2].split(',')
        kmers = set()
        for h in haps:
            if h == 'dm': continue
            c = h.split('_')[0]
            i = int(h.split('hap')[1])
            try:
                with open(f'{cwd}/kmer_size_{k}/{c}_hap{i-4}/syn_fasta_{c}_hap{i}_bin_{b}.k{k}.good.txt', 'r') as f:
                    hkmers = set([l.strip().split()[0] for l in f])
                kmers.update(hkmers)
            except FileNotFoundError:
                pass
        return kmers
    # END

    def getwinmers(grp):
        kmers = set()
        for row in grp[1].itertuples(index=False):
            kmers.update(gethapmers(row))
        return kmers

    hapdf = pd.read_table('/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/haplotype_graph.txt', header=None)


    # <editor-fold desc="Find kmers that are unique in each node: Get kmers that are present in a single node in the graph">
    unikmers = set()
    badkmers = set()
    i = 1
    for row in tqdm(hapdf.itertuples(index=False)):
        kmers = gethapmers(row)
        kmers.difference_update(badkmers)
        badkmers.update(kmers.intersection(unikmers))
        unikmers.symmetric_difference_update(kmers)

    with open(f'{cwd}/kmer_size_{k}/uninodekmers.txt', 'w') as fout:
        fout.write('\n'.join(unikmers))

    # Get kmers that are unique in each node and save them
    unikmers = set([l.strip() for l in open(f'{cwd}/kmer_size_{k}/uninodekmers.txt')])
    with open(f'{cwd}/kmer_size_{k}/nodekmers.txt', 'w') as fout:
        for row in tqdm(hapdf.itertuples(index=False)):
            kmers = gethapmers(row)
            kmers.intersection_update(unikmers)
            unikmers.difference_update(kmers)
            if len(kmers) > 0:
                fout.write("\n".join([f'{row[3]}\t{k}' for k in kmers]) + "\n")
    # </editor-fold>

    # <editor-fold desc="Find kmers that are unique in to each window (but might be shared between nodes): Get kmers that are present in a single window in the graph">

    unikmers = set()
    badkmers = set()
    i = 1
    for grp in tqdm(hapdf.groupby(0)):
        kmers = getwinmers(grp)
        kmers.difference_update(badkmers)
        badkmers.update(kmers.intersection(unikmers))
        unikmers.symmetric_difference_update(kmers)

    with open(f'{cwd}/kmer_size_{k}/uniwinkmers.txt', 'w') as fout:
        fout.write('\n'.join(unikmers))

    # Get kmers that are unique in each node and save them
    unikmers = set([l.strip() for l in open(f'{cwd}/kmer_size_{k}/uniwinkmers.txt')])
    with open(f'{cwd}/kmer_size_{k}/winkmers.txt', 'w') as fout:
        for grp in tqdm(hapdf.groupby(0)):
            winkmers = getwinmers(grp)
            for row in grp[1].itertuples(index=False):
                kmers = gethapmers(row)
                kmers.intersection_update(unikmers)
                if len(kmers) > 0:
                    fout.write("\n".join([f'{row[3]}\t{k}' for k in kmers]) + "\n")
            unikmers.difference_update(winkmers)
    # </editor-fold>

    return
# END


def summary_plots_kmers_per_node():
    """
    Collection of plots to describe how many nodes have kmers and various distributions
    :return:
    """
    from collections import deque, Counter
    from matplotlib import pyplot as plt
    from tqdm import tqdm
    from hometools.plot import cleanax

    allnodes = tuple(range(2927))
    nodecnts = deque()
    i = 100
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


def count_nodekmers_from_samples(samplekmersfin):
    """
    From kmer set for a sample, select node-specific kmers.
    :return:
    """

    samplekmersfin = '/home/ra98jam/projects/pantohap/results/potato_haps/otava_k51_filter.good.kmers.gz'
    fo = gzopen if isgzip(samplekmersfin) else open

    with open('/home/ra98jam/projects/pantohap/results/potato_haps/nodekmers.txt', 'r') as fin:
        kmers = {line.strip().split()[1].encode(): 0 for line in fin}
    # Get kmer counts for sample
    with fo(samplekmersfin, 'r') as fin:
        for line in tqdm(fin):
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
    with open('/home/ra98jam/projects/pantohap/results/potato_haps/nodekmers.txt', 'r') as fin:
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
    with open('/home/ra98jam/projects/pantohap/results/potato_haps/node_k_stats.txt', 'w') as fout:
        for k in nodekmercnts.keys():
            m = round(np.mean(nodekmercnts[k]), 2) if len(nodekmercnts[k]) > 0 else 0
            fout.write(f'{k}\t{nodeks[k]}\t{len(nodekmercnts[k])}\t{round(len(nodekmercnts[k])/nodeks[k], 2)}\t{m}\t{",".join(list(map(str, nodekmercnts[k])))}\n')

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






    nodestats = nodestats.loc[nodestats[2] >= 0.8]


    # Cound number of kmers for each

    return
# END

