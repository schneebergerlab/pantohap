# <editor-fold desc='Step1: Generate haplotype graph'>
# Update SNPs genotype to include the information about the deletion or SR markers
# SNP_type - Genotype
# REF_allele in syntenic            -           0
# ALT_allele in syntenic            -           1
# Deleted                           -           2
# In inversion                      -           3
# In TDs (considered as deletion)   -           2


# Currently use the graph generated using binning of reference genomes
util.hapnodesfromvcffile()

# TODO: Consider using msyd hap-graph for this problem
run_msyd.sh
# </editor-fold>

# <editor-fold desc='Kmer counting'>
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl

# For each genome, select kmers present only once in the assembly
get_unique_kmers_from_assemblies_part1.sh

# Merge multimers
# SH
for k in 21 31 41 51; do
    cd ${cwd}/kmer_size_${k}
    $meryl union output all_genome_multi */*multi &
done

get_unique_kmers_from_assemblies_part2.sh

# Merge unimers
# SH
for k in 21 31 41 51; do
    cd ${cwd}/kmer_size_${k}
    $meryl union output all_genome_good */*good &
done


# for each 100kb window, get syntenic sequence from query (40 hap) genomes.
get_node_query_sequence()

# Get kmers for each of the fasta file created using get_node_query_sequence()
SBATCH_get_unique_kmer_per_window.sh

# Merge Kmers from each collapsed node and then find marker kmers (that are present in that node only) for each node
from functools import partial
from multiprocessing import Pool
cwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/'
cwd = '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/'
with Pool(processes=4) as pool:
    pool.map(partial(get_unique_kmers_per_node, cwd=cwd), [21, 31, 41, 51])


# Generate summary statistics
summary_plots_kmers_per_node()

# Reads kmers from a sample of interest and get kmer counts for the nodekmers


# </editor-fold>
