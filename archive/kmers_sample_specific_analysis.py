import socket
from string import ascii_uppercase
import pysam
from multiprocessing import Pool
from collections import defaultdict
import pandas as pd
import seaborn as sns

GENOMES = list(ascii_uppercase[:10])
HAPS = list(range(1, 5))

if socket.gethostname() == 'LMBIDP1-LXSCH01':
    # Indir for local compute
    CWD = '/home/ra98jam/d16/projects/potato_hap_example//results/kmer_analysis/marker_kmers_21/'
else:
    # Indir for cluster compute
    CWD = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example//results/kmer_analysis/marker_kmers_21/'

def get_chr_counts(gen):
    hap_cnts = dict()
    for hap in HAPS:
        chr_counts = defaultdict(int)
        with pysam.AlignmentFile(f'{CWD}/{gen}_hap{hap}/{gen}_hap{hap}_uni.sorted.bam') as bam:
            for read in bam:
                chr_counts[read.reference_name] += 1
        keys = list(chr_counts.keys())
        for k in keys:
            if "scaffold" in k:
                chr_counts.pop(k)
        hap_cnts[hap] = {k.split('_')[0]: v for k,v in chr_counts.items()}
    return hap_cnts
# END

with Pool(processes=10) as pool:
    gen_chr_counts = pool.map(get_chr_counts, GENOMES)
gen_chr_counts_df = dict(zip(GENOMES, gen_chr_counts))
gen_chr_counts_df = pd.concat({k: pd.DataFrame(v) for k, v in gen_chr_counts_df.items()})
gen_chr_counts_df.reset_index(drop=False, inplace=True)
gen_chr_counts_df.columns = 'genome chromosome hap1 hap2 hap3 hap4'.split()
gen_chr_counts_df = gen_chr_counts_df.melt(id_vars='genome chromosome'.split())
gen_chr_counts_df.columns = 'genome chromosome haplotype counts'.split()
gen_chr_counts_df.plot(kind='scatter', x='chromosome', y='counts')

ax = sns.stripplot(data=gen_chr_counts_df, x='chromosome', y='counts', hue='genome')
ax.set_yscale('log')
ax = sns.violinplot(data=gen_chr_counts_df, x='chromosome', y='counts', alpha=0.5, ax=ax)
plt.tight_layout()
plt.legend(ncol=5)
plt.savefig(f'{CWD}kmers_sample_specific_counts.pdf')




# <editor-fold desc="SCRATCH AND TESTS>
gc = deque()
with open("/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/node_kmers/chr01/nodekmers_k51_2024_08_16_haplotype_graph_chr01_div0.1_2024_08_02.txt") as f:
    for line in tqdm(f):
        kmer = line.strip().split()[1]
        gc.append((kmer.count('G') + kmer.count('C'))/len(kmer))
        break

# </editor-fold>