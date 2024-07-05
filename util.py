# <editor-fold desc='Define imports'>
import gzip
import string
from collections import defaultdict, deque
from itertools import product
from functools import partial
from multiprocessing import Pool
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import igraph as ig
from gzip import open as gzopen
from tqdm import tqdm
import seaborn as sns

from hometools.plot import cleanax, loghist
from hometools.hometools import isgzip, mapbp

# </editor-fold>


# <editor-fold desc='Define constants'>
chrsize = 46102915
samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]        # Without Otava
samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10]+'O', range(5, 9))]    # With Otava
# </editor-fold>


# <editor-fold desc='Classes'>

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
    return sum([i != j for i, j in zip(iter1, iter2)])
# END


def getsamplestates(sid, pwd, snpdf):
    print(sid)
    # TODO: Update input file name
    hdf = pd.read_table(f'{pwd}/dm_{sid[0]}_chr02_hap{sid[1]}syri.out.bed.snp_anno.txt', header=None, low_memory=False)
    # hdf = pd.read_table(f'{pwd}/dm_{sid[0]}_chr02_hap{sid[1]}.with_Otava.syri.out.bed.snp_anno.txt', header=None, low_memory=False)
    hdf = hdf.loc[hdf[6] != 'SNP']
    hdf.sort_values([1, 2], inplace=True)
    state = deque()
    colindex = [i for i, c in enumerate(snpdf.columns) if c == sid][0]
    for row in tqdm(snpdf.itertuples()):
        # Condition 0: Syntenic SNP already found
        if row[colindex + 1] == 1:
            state.append(1)
            continue
        hdfsnp = hdf.loc[(hdf[1] <= int(row[0][1])) & (hdf[2] >= int(row[0][1]))]
        # Condition 1: Is syntenic but no SNP was found
        if 'SYNAL' in hdfsnp[6].values:
            state.append(0)
            continue
        # Condition 2: Has an inverted allele
        if 'INVAL' in hdfsnp[6].values:
            state.append(2)
            continue
        # Condition 4: Everything else is considered as a deletion
        state.append(3)
    return state
# END


def updatesnpgenotypes():
    from tqdm import tqdm
    from itertools import product

    """
        # SNP_type - Genotype
        # REF_allele in syntenic            -           0
        # ALT_allele in syntenic            -           1
        # In inversion                      -           2
        # Deleted                           -           3
        # In TDs (considered as deletion)   -           3
    """

    # Get the syri annotations at the identified SNP marker positions
    """ SH
    # Without Otava
    fname=dm_all_sample_chr2.syri.nosr.snps.merged.vcf.regions_for_tabix.txt
    awk 'NR>1 { print $1"\t"$2-1"\t"$2+1}' dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt > $fname
    for c in {A..J}; do
        for h in {5..8}; do
            tabix dm_${c}_chr02_hap${h}syri.out.bed.gz -R $fname | sort | uniq > dm_${c}_chr02_hap${h}syri.out.bed.snp_anno.txt &
        done
        wait
    done
    
    # With Otava
    fname=dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.regions_for_tabix.txt
    awk 'NR>1 { print $1"\t"$2-1"\t"$2+1}' dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.txt > $fname
    for c in {A..J} O; do
        {
        for h in {5..8}; do
            tabix dm_${c}_chr02_hap${h}syri.out.bed.gz -R $fname | sort | uniq > dm_${c}_chr02_hap${h}.with_Otava.syri.out.bed.snp_anno.txt
        done
        } &
    done
    """

    pwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data'
    # TODO: Update the input file name
    snpfin = f"dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt"
    # snpfin = f"dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.txt"

    sids = list(product(string.ascii_uppercase[:10], range(5, 9)))
    # sids = list(product(string.ascii_uppercase[:10]+'O', range(5, 9)))
    # Read snpfin
    snpdf = dict()
    with open(f'{pwd}/{snpfin}', 'r') as fin:
        _ = fin.readline()
        for i, line in tqdm(enumerate(fin)):
            # if i == 10: break
            line = line.strip().split()
            v = np.array([int(i) if i != '.' else 0 for i in line[7:]])
            if set(v) == set([0]): continue
            snpdf[tuple(line[:4])] = v
    snpdf = pd.DataFrame.from_dict(snpdf, orient='index')
    snpdf.columns = sids

    with Pool(processes=20) as pool:
        states = pool.map(partial(getsamplestates, pwd=pwd, snpdf=snpdf), sids)

    snpdf.columns = list(range(len(sids)))
    for i, state in enumerate(states):
        snpdf[i] = state

    anndf = pd.DataFrame([list(i) for i in snpdf.index.values])
    anndf.columns = '#CHROM POS REF ALT'.split()
    anndf['QUAL'] = '.'
    anndf['FILTER'] = 'PASS'
    anndf['FORMAT'] = 'GT'
    snpdf.columns = ["_hap".join(list(map(str, sid))) for sid in sids]
    snpdfout = snpdf.reset_index(drop=True)
    snpdfout = pd.concat([anndf, snpdfout], axis=1)

    # TODO: Update the output file name
    snpdfout.to_csv(f'{pwd}/dm_all_sample_chr2.syri.nosr.snps.merged.vcf.del_markers.txt', sep='\t', index=False)
    # snpdfout.to_csv(f'{pwd}/dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.del_markers.txt', sep='\t', index=False)

    return
# END


def hapnodesfromvcffile(snpfin):
    import pandas as pd
    import numpy as np
    from collections import deque
    from itertools import product
    from tqdm import tqdm
    import igraph as ig
    from datetime import datetime
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

    def getnodes(start):
        print(start)
        # end = a[i + 1]
        end = start + 100000    # Window size
        snps = snpdf.loc[(snpdf.position >= start) & (snpdf.position < end)].copy()
        # snpdf.drop(snps.index, inplace=True)
        snpgt = snps.iloc[:, 7:]
        snpgt = snpgt.T
        if snpgt.shape[1] == 0:
            return [start, end]

        haplist = [grp for i, grp in snpgt.groupby(list(snpgt.columns))]
        nhap = len(haplist)

        # Find groups of haps that are very similar (divergence < 10%)
        # dist = [editdist(i.iloc[0], j.iloc[0]) for i, j in product(haplist, haplist)]
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
        hapsample = [','.join(sorted(h.index.values)) for h in hapglist]
        return [start, end] + hapsample
    # END

    # snpfin = f"/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt"
    pwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data'
    pwd = '/home/ra98jam/d16/projects/potato_hap_example/data/'
    snpfin = 'dm_all_sample_chr2.syri.nosr.snps.merged.vcf.del_markers.txt'
    snpfin = 'dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.del_markers.txt'

    snpdf = pd.read_table(f"{pwd}/{snpfin}")
    snpdf.columns = 'chr position reference alternate'.split() + list(snpdf.columns[4:])
    snpdf['dm'] = '.'
    snpdf.replace('.', 0, inplace=True)
    snpdf.replace('1', 1, inplace=True)
    snpdf.replace('2', 2, inplace=True)
    snpdf.replace('3', 3, inplace=True)

    a = np.arange(1, chrsize, 100000)
    a = np.append(a, chrsize)

    with Pool(processes=20) as pool:
        hapblocks = pool.map(getnodes, a[:-1])

    # hapblocks = deque()
    # for i, start in tqdm(enumerate(a[:-1])):
    #     end = a[i + 1]
    #     snps = snpdf.loc[(snpdf.position >= start) & (snpdf.position < end)].copy()
    #     snpdf.drop(snps.index, inplace=True)
    #     snpgt = snps.iloc[:, 7:]
    #     snpgt = snpgt.T
    #     if snpgt.shape[1] == 0:
    #         continue
    #
    #     haplist = [grp for i, grp in snpgt.groupby(list(snpgt.columns))]
    #     nhap = len(haplist)
    #
    #     # Find groups of haps that are very similar (divergence < 10%)
    #     dist = [editdist(i.iloc[0], j.iloc[0]) for i, j in product(haplist, haplist)]
    #     dmat = np.reshape([editdist(i.iloc[0], j.iloc[0]) for i, j in product(haplist, haplist)], [nhap, nhap])
    #     conn = (dmat / haplist[0].shape[1]) < 0.10
    #     conn = np.triu(conn, 1)
    #     g = ig.Graph.Adjacency(conn)
    #     hapgrp = g.connected_components(mode='weak')
    #
    #     # Merge haps that are within the same group
    #     hapglist = []
    #     for hg in hapgrp:
    #         hap = pd.concat([haplist[i] for i in hg])
    #         hapglist.append(hap)
    #     nhap = len(hapglist)
    #     hapsample = [','.join(sorted(h.index.values)) for h in hapglist]
    #     hapblocks.append([start, end] + hapsample)

    # hapoblist = deque()
    # cnt = 0
    # for i, h in enumerate(hapblocks):
    #     for ss in h[2:]:
    #         ss = sorted(ss.split(','))
    #         hapoblist.append(hapobject(cnt, h[0], h[1] - 1, ss))
    #         cnt += 1
    index = 0
    # with open(f"{pwd}/haplotype_graph_{str(datetime.now().date()).replace('-', '_')}.txt", 'w') as fout:
    with open(f"{pwd}/haplotype_graph_with_Otava{str(datetime.now().date()).replace('-', '_')}.txt", 'w') as fout:
        for hapblock in hapblocks:
            for hap in hapblock[2:]:
                fout.write(f'{hapblock[0]}\t{hapblock[1]}\t{hap}\t{index}\n')
                index += 1

    return
# END


def get_sequence(sample, pwd, a):
        print(sample)
        g, h = sample.split('_hap')
        # get syn regions
        synr = pd.read_csv(f'{pwd}/dm_{g}_chr02_hap{h}syri.out', header=None, sep='\t')
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
        (chrom, seq) = next(readfasta_iterator(open(f'{pwd}/{g}_chr02_hap{h}.fa', 'r')))
        TOBREAK = False
        for i, r in synd.items():
            with open(f'{pwd}/synfastas/syn_fasta_{sample}_bin_{i}.fasta', 'w') as fout:
                for c in r.df.itertuples(index=False):
                    # Get query start position
                    if c.Start != i:
                        qs = c[3][0]
                    else:
                        # Get the exact matching position when the reference alignment was (potentially) trimmed at the start of the window
                        mappos = mapbp(sfin=f'{pwd}/dm_{g}_chr02_hap{h}syri.out.bed.gz', mapfin=f'{pwd}/dm_{g}_chr02_hap{h}.bam', d=True, posstr=f'chr02:{i}-{i}')
                        mappos = list(map(lambda x: int(x.split(':')[1].split('-')[0]) if '+' in x else None, mappos))
                        mappos = [x for x in mappos if x is not None]
                        # If more than 1 mapping position was identified, then select the position that is within the selected alignment
                        if len(mappos) > 1:
                            mapfit = list(map(lambda x: c[3][0] < x < c[3][1], mappos))
                            if mapfit.count(True) == 1:
                                qs = mappos[mapfit.index(True)]
                            else:
                                print('Breaking qs')
                                TOBREAK = True
                                break
                        else:
                            qs = mappos[0]

                    # Get query end position
                    if c.End != i+100000:
                        qe = c[3][1]
                    else:
                        # Get the exact matching position when the reference alignment was (potentially) trimmed at the end of the window
                        mappos = mapbp(sfin=f'{pwd}/dm_{g}_chr02_hap{h}syri.out.bed.gz', mapfin=f'{pwd}/dm_{g}_chr02_hap{h}.bam', d=True, posstr=f'chr02:{i+100000}-{i+100000}')
                        mappos = list(map(lambda x: int(x.split(':')[1].split('-')[0]) if '+' in x else None, mappos))
                        mappos = [x for x in mappos if x is not None]
                        # If more than 1 mapping position was identified, then select the position that is within the selected alignment
                        if len(mappos) > 1:
                            mapfit = list(map(lambda x: c[3][0] < x < c[3][1], mappos))
                            if mapfit.count(True) == 1:
                                qe = mappos[mapfit.index(True)]
                                # print(c, mappos, qe)
                            else:
                                print('Breaking qe')
                                TOBREAK = True
                                break
                        else:
                            qe = mappos[0]
                    fout.write(f'>chr02_{g}_hap{h}_r_{c[1]}_{c[2]}_q_{qs}_{qe}\n')
                    fout.write(seq[(qs-1):qe] + '\n')
                if TOBREAK:
                    break


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
    from multiprocessing import Pool
    from functools import partial



    chrsize = 46102915
    samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]
    pwd = '/home/ra98jam/d16/projects/potato_hap_example/data/'
    a = np.arange(1, chrsize, 100000)
    a = np.append(a, chrsize)
    with Pool(processes=10) as pool:
        pool.map(partial(get_sequence, pwd=pwd, a=a), samples[1:])
    return

    #
    # # Run in /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/
    # # For each sample, select syntenic region in each window (100kb). Skip DM.
    # for sample in samples[1:5]:
    #     print(sample)
    #     g, h = sample.split('_hap')
    #     # get syn regions
    #     synr = pd.read_csv(f'{pwd}/dm_{g}_chr02_hap{h}syri.out', header=None, sep='\t')
    #     synr = synr.loc[synr[10] == 'SYNAL']
    #     synr[[1, 2, 6, 7]] = synr[[1, 2, 6, 7]].astype(int)
    #     synd = dict()
    #     for i in a:
    #         sinr = synr.loc[(synr[1] < (i+100000)) & (synr[2] > i)]
    #         if sinr.shape[0] > 0:
    #             pr1 = pr.from_dict({"Chromosome": ["chr02"],
    #                                 "Start": [i],
    #                                 "End": [i+100000]})
    #             pr2 = pr.from_dict({"Chromosome": ["chr02"]*sinr.shape[0],
    #                                 "Start": sinr[1],
    #                                 "End": sinr[2],
    #                                 "query": sinr[[6, 7]].values.tolist()})
    #             synd[i] = pr2.intersect(pr1)
    #     (chrom, seq) = next(readfasta_iterator(open(f'{pwd}/{g}_chr02_hap{h}.fa', 'r')))
    #
    #     TOBREAK = False
    #     for i, r in synd.items():
    #         with open(f'{pwd}/synfastas/syn_fasta_{sample}_bin_{i}.fasta', 'w') as fout:
    #             for c in r.df.itertuples(index=False):
    #                 # Get query start position
    #                 if c.Start != i:
    #                     qs = c[3][0]
    #                 else:
    #                     # Get the exact matching position when the reference alignment was (potentially) trimmed at the start of the window
    #                     mappos = mapbp(sfin=f'{pwd}/dm_{g}_chr02_hap{h}syri.out.bed.gz', mapfin=f'{pwd}/dm_{g}_chr02_hap{h}.bam', d=True, posstr=f'chr02:{i}-{i}')
    #                     mappos = list(map(lambda x: int(x.split(':')[1].split('-')[0]) if '+' in x else None, mappos))
    #                     mappos = [x for x in mappos if x is not None]
    #                     # If more than 1 mapping position was identified, then select the position that is within the selected alignment
    #                     if len(mappos) > 1:
    #                         mapfit = list(map(lambda x: c[3][0] < x < c[3][1], mappos))
    #                         if mapfit.count(True) == 1:
    #                             qs = mappos[mapfit.index(True)]
    #                         else:
    #                             print('Breaking qs')
    #                             TOBREAK = True
    #                             break
    #                     else:
    #                         qs = mappos[0]
    #
    #                 # Get query end position
    #                 if c.End != i+100000:
    #                     qe = c[3][1]
    #                 else:
    #                     # Get the exact matching position when the reference alignment was (potentially) trimmed at the end of the window
    #                     mappos = mapbp(sfin=f'{pwd}/dm_{g}_chr02_hap{h}syri.out.bed.gz', mapfin=f'{pwd}/dm_{g}_chr02_hap{h}.bam', d=True, posstr=f'chr02:{i+100000}-{i+100000}')
    #                     mappos = list(map(lambda x: int(x.split(':')[1].split('-')[0]) if '+' in x else None, mappos))
    #                     mappos = [x for x in mappos if x is not None]
    #                     # If more than 1 mapping position was identified, then select the position that is within the selected alignment
    #                     if len(mappos) > 1:
    #                         mapfit = list(map(lambda x: c[3][0] < x < c[3][1], mappos))
    #                         if mapfit.count(True) == 1:
    #                             qe = mappos[mapfit.index(True)]
    #                             # print(c, mappos, qe)
    #                         else:
    #                             print('Breaking qe')
    #                             TOBREAK = True
    #                             break
    #                     else:
    #                         qe = mappos[0]
    #                 fout.write(f'>chr02_{g}_hap{h}_r_{c[1]}_{c[2]}_q_{qs}_{qe}\n')
    #                 fout.write(seq[(qs-1):qe] + '\n')
    #             if TOBREAK:
    #                 break
    return
# END


def get_unique_kmers_per_node(k, cwd):
    """
    Read the haplotype graph and:
    1) Get set of kmers for each node
    2) Filter out kmers that are present in more than 1 node
    3) Save the unique kmers for each node
    :return:
    """
    from datetime import datetime
    getdate = lambda: str(datetime.now().date()).replace('-', '_')

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

    # hapdf = pd.read_table('/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/haplotype_graph.txt', header=None)
    hapdf = pd.read_table(f'{cwd}../../data/haplotype_graph_2024_06_14.txt', header=None)

    # <editor-fold desc="Find kmers that are unique in each node: Get kmers that are present in a single node in the graph">
    unikmers = set()
    badkmers = set()
    i = 1
    for row in tqdm(hapdf.itertuples(index=False)):
        kmers = gethapmers(row)
        kmers.difference_update(badkmers)
        badkmers.update(kmers.intersection(unikmers))
        unikmers.symmetric_difference_update(kmers)

    with open(f'{cwd}/kmer_size_{k}/uninodekmers_{getdate()}.txt', 'w') as fout:
        fout.write('\n'.join(unikmers))

    # Get kmers that are unique in each node and save them
    unikmers = set([l.strip() for l in open(f'{cwd}/kmer_size_{k}/uninodekmers_{getdate()}.txt', 'r')])
    with open(f'{cwd}/kmer_size_{k}/nodekmers_{getdate()}.txt', 'w') as fout:
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

    with open(f'{cwd}/kmer_size_{k}/uniwinkmers_{getdate()}.txt', 'w') as fout:
        fout.write('\n'.join(unikmers))

    # Get kmers that are unique in each node and save them
    unikmers = set([l.strip() for l in open(f'{cwd}/kmer_size_{k}/uniwinkmers_{getdate()}.txt')])
    with open(f'{cwd}/kmer_size_{k}/winkmers_{getdate()}.txt', 'w') as fout:
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


def get_threads():
    import igrap

    def find_all_paths(graph, start, end):
        def dfs(current, path):
            path.append(current)
            if current == end:
                all_paths.append(list(path))
            else:
                for neighbor in graph[current]:
                    if neighbor not in path:
                        dfs(neighbor, path)
            path.pop()

        all_paths = []
        dfs(start, [])
        return all_paths

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

    cwd = '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/'
    emoutfin = "/home/ra98jam/d16/test/WhiteRose_results_2024_06_14/EM_results/EM.v03.WhR.w0.results.tsv"
    hapfin = "../data/haplotype_graph_2024_06_14.txt"

    hapoblist = deque()
    with open(hapfin, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            hapoblist.append(hapobject(int(line[3]), int(line[0]), int(line[1]), sorted(line[2].split(','))))

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

    # Add ploidy levels in the graph
    ploidy = deque()
    emout = pd.read_csv(emoutfin, header=None, sep='\t')
    def get_ploidy(row):
        print(row[9:])
        if all([r == 0 for r in row[9:]]):
            return 0
        # TODO: for a node, I assign based on the maximum probablity from EM.
        #  So, if for a node, EM outputs probablity of 51% for 1 copy and 49%
        #  for 2 copy, I assign ploidy of 1 to the node. Probably, this can be
        #  fine-tuned.
        return row[9:].to_list().index(max(row[9:]))
    # END

    ploidy = emout.apply(get_ploidy, axis=1)
    G.vs['ploidy'] = ploidy

    # Delete edges between nodes with unequal ploidy
    to_del = [e for e in G.get_edgelist() if G.vs[e[0]]['ploidy'] != G.vs[e[1]]['ploidy']]
    G.delete_edges(to_del)

    # Delete vertices with ploidy 0
    G.delete_vertices([i for i, p in enumerate(ploidy) if p == 0])

    # Check whether the predicted ploidy is correct
    ploidy_color = {
        0: "lightgrey",
        1: "orange",
        2: "blue",
        3: "red",
        4: "black"
    }
    cor_plo = deque()   # Green for correct ploidy, white for wrong
    for v in G.vs:
       p = v['ploidy']
       gs = vars(hapoblist[v['name']])['genomes']
       if p == len([g for g in gs if 'A_' in g ]):
           cor_plo.append(ploidy_color[v['ploidy']])
       else:
           cor_plo.append(ploidy_color[0])
    G.vs['cor'] = cor_plo

    # Save figures
    ig.plot(G, target='tmp.pdf', layout='kk', vertex_size=np.array(G.degree())+4, vertex_color=G.vs['cor'], vertex_frame_color="black", vertex_frame_width=0,  edge_width=0.1, edge_arrow_size=0.2, bbox=[800, 800])
# END

# <editor-fold desc="OBSOLETE FUNCTIONS">

def count_kmers_from_samples(samplekmersfin):
    """
    For a given sample, fetch kmer counts for node-specific and window-specific kmers.
    :return: Generates
    """
    import os
    from hometools.hometools import isgzip
    from gzip import open as gzopen
    from collections import defaultdict, deque
    import numpy as np
    import pandas as pd
    from tqdm import tqdm

    ks = [21, 31, 41, 51]
    cwd = '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/'

    for k in ks:
        os.chdir(f'{cwd}/kmer_size_{k}/')
        samplekmersfin = f'./otava_kmers/otava_k{k}_filter.good.kmers.gz'
        fo = gzopen if isgzip(samplekmersfin) else open

        with open('nodekmers.txt', 'r') as fin:
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
        with open('nodekmers.txt', 'r') as fin:
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
        with open('node_k_stats.txt', 'w') as fout:
            for k in nodekmercnts.keys():
                m = round(np.mean(nodekmercnts[k]), 2) if len(nodekmercnts[k]) > 0 else 0
                fout.write(f'{k}\t{nodeks[k]}\t{len(nodekmercnts[k])}\t{round(len(nodekmercnts[k])/nodeks[k], 2)}\t{m}\t{",".join(list(map(str, nodekmercnts[k])))}\n')

        nodestats = pd.read_table('node_k_stats.txt', header=None)
        nodestats.sort_values([2, 3], ascending=False, inplace=True)
        nodestats.to_csv('node_k_stats.txt', sep='\t', header=False, index=False)

    # Generate plots here: analysis_plots.node_k_stats_plots
    return
# END


# </editor-fold>

