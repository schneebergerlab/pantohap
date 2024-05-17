# <editor-fold desc='Step1: Generate haplotype graph'>

# Currently use the graph generated using binning of reference genomes

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
with Pool(processes=4) as pool:
    for k in [21, 31, 41, 51]:
        get_unique_kmers_per_node(cwd, k) # Creates nodekmers.txt


# Generate summary statistics
summary_plots_kmers_per_node()

# Reads kmers from a sample of interest and get kmer counts for the nodekmers


# </editor-fold>
