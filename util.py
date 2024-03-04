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


    snpfin = f"/home/ra98jam/projects/pantohap/data/dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt"
    snpfin = f"dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt"
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




        # TODO: Use the syntenic regions in synd to fetch sequence for syntenic region in each hapblock

    for i, h in enumerate(hapblocks):
        rs = h[0]
        re = h[1]
        for j, ss in enumerate(h[2:]):
            # Open fasta file where all the fasta sequences for a hapblock would be saves
            with open(f'hapblocksequence/hap{i}_block{j}.fa', 'r') as fout:
                for s in ss:
                    if s == 'dm':






    hapoblist = deque()
    cnt = 0
    for i, h in enumerate(hapblocks):
        for ss in h[2:]:
            ss = sorted(ss.split(','))
            hapoblist.append(hapobject(cnt, h[0], h[1] - 1, ss))
            cnt += 1
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