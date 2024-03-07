from collections import defaultdict

import numpy as np
import pandas as pd


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

def get_unique_kmers_per_node():
    """
    Read the haplotype graph and:
    1) Get set of kmers for each node
    2) Filter out kmers that are present in more than 1 node
    3) Save the unique kmers for each node
    :return:
    """

    from tqdm import tqdm
    from hometools.hometools import view
    import pandas as pd
    hapdf = pd.read_table('/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/haplotype_graph.txt', header=None)

    # Get kmers that are present in a single node in the graph
    unikmers = set()
    badkmers = set()
    i = 1
    for row in tqdm(hapdf.itertuples(index=False)):
        s = row[0]
        haps = row[2].split(',')
        kmers = set()
        for h in haps:
            try:
                hkmers = set(pd.read_table(f'synfastas/syn_fasta_{h}_bin_{s}.k51.unique.txt', engine='c', header=None)[0])
                kmers.update(hkmers)
            except FileNotFoundError:
                pass
        kmers.difference_update(badkmers)
        badkmers.update(kmers.intersection(unikmers))
        unikmers.symmetric_difference_update(kmers)
    with open('/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/unikmers.txt', 'w') as fout:
        fout.write('\n'.join(unikmers))


    # Get kmers that are unique in each node and save them
    unikmers = set([l.strip() for l in open('unikmers.txt')])
    with open('nodekmers.txt', 'w') as fout:
        for row in tqdm(hapdf.itertuples(index=False)):
            s = row[0]
            haps = row[2].split(',')
            kmers = set()
            for h in haps:
                try:
                    hkmers = set(pd.read_table(f'synfastas/syn_fasta_{h}_bin_{s}.k51.unique.txt', engine='c', header=None)[0])
                    kmers.update(hkmers)
                except FileNotFoundError:
                    pass
            kmers.intersection_update(unikmers)
            unikmers.difference_update(kmers)
            fout.write("\n".join([f'{row[3]}\t{k}' for k in kmers]) + "\n")
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
    plt.close()



