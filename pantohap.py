# <editor-fold desc='Step1: Generate haplotype graph'>
# Update SNPs genotype to include the information about the deletion or SR markers
# SNP_type - Genotype
# REF_allele in syntenic            -           0
# ALT_allele in syntenic            -           1
# Deleted                           -           2
# In inversion                      -           3
# In TDs (considered as deletion)   -           2


# Generate haplotype genotype file required for generating haplotype graph
get_haplotype_vcf.sh

# Currently the graph is generated using binning of reference genomes (util.py)
util.hapnodesfromvcffile()

# TODO: Consider using msyd hap-graph for this problem
run_msyd.sh
# </editor-fold>

# <editor-fold desc='Kmer counting'>
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl

# For each genome, select kmers present only once in the assembly
sbatch ../../tool/sbatch_files/SBATCH_get_unique_kmers_from_assemblies_part1.sh


# Merge multimers
# SH
# for k in 21 31 41 51; do
for k in 51; do
    cd ${cwd}/kmer_size_${k}
    $meryl union output all_genome_multi *hap*/*multi
done

sbatch ../../tool/sbatch_files/SBATCH_get_unique_kmers_from_assemblies_part2.sh

# Merge unimers
# SH
# for k in 21 31 41 51; do
for k in 51; do
    cd ${cwd}/kmer_size_${k}
    $meryl union output all_genome_good *hap*/*good &
done


# for each 100kb window, get syntenic sequence from query (40 hap) genomes.
cwd=/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/node_kmers/
util.get_node_query_sequence()

# Get kmers for each of the fasta file created using get_node_query_sequence()
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/node_kmers/
for c in {01..12}; do
    cd $cwd ; cd chr$c
    chrid=$c
    export chrid
    sbatch ../../../../tool/sbatch_files/SBATCH_get_unique_kmer_per_window.sh
done

# Merge Kmers from each collapsed node and then find marker kmers (that are present in that node only) for each node
from functools import partial
from multiprocessing import Pool
cwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/'
cwd = '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/'
with Pool(processes=4) as pool:
    pool.map(partial(get_unique_kmers_per_node, cwd=cwd), [21, 31, 41, 51])

# Haplotype graph filename
k=51
chrs = ["chr{:02d}".format(i) for i in range(1, 13)]
for d in [0.1, 0.05, 0.01]:
for d in [0.05, 0.01]:
    hgf = 'haplotype_graph_{c}' + '_div{d}_2024_08_02.txt'.format(d=d)
    with Pool(processes=12) as pool:
        pool.map(partial(get_unique_kmers_per_node, hgf=hgf, k=k), chrs)


# Run EM to get ploidy levels for each node

# Get a connected paths through the graph



# Generate summary statistics
summary_plots_kmers_per_node()

# Reads kmers from a sample of interest and get kmer counts for the nodekmers


# </editor-fold>
