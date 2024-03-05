# <editor-fold desc='Step1: Generate haplotype graph'>

# Currently use the graph generated using binning of reference genomes

# </editor-fold>


# <editor-fold desc='Kmer counting'>

# for each 100kb window, get syntenic sequence from query (40 hap) genomes.
get_node_query_sequence()

# Get kmers for each of the fasta file created using get_node_query_sequence()
SBATCH_get_unique_kmer_per_window.sh

# Merge Kmers from each collapsed node and then find marker kmers (that are present in that node only) for each node
get_unique_kmers_per_node()

# Generate summary statistics


# </editor-fold>
