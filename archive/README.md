# README
Pantohap generates a haplotype graph by finding shared haplotypes in multiple genomes.

The steps of the pipeline are present in the pantohap.sh. The steps are:

1) Compare the genomes to the reference and identify shared variations between them. 
   1) Done using the get_haplotype_vcf.sh and involves:
      1) Aligning genomes to reference
      2) Variant calling
      3) Get short variants in syntenic regions
      4) Merge and normalise variant calls based on structural properties (e.g. large deletions, genomic rearrangements, etc)
2) Use the variants identified above to create haplotype graph. For this, the reference genome is divided into 100kbp windows and genomes sharing similar variants in a window are collapsed to form nodes. Done using the util.py -> hapnodesfromvcffile() function. 
3) Get candidate marker kmers from the genome assemblies. Kmers are identified in the genome assemblies and any kmer that is present more than once in any genome is filtered out.
   1) Uses sbatch_files/SBATCH_get_unique_kmers_from_assemblies_part1.sh and sbatch_files/SBATCH_get_unique_kmers_from_assemblies_part2.sh files. 
4) For each window, fetch the fasta sequence of the corresponding syntenic region for each of the genomes. Done using util.py -> get_node_query_sequence()
5) From these fasta files, select kmers that are present only once in fasta and selected as candidate kmers from the assemblies. Done using sbatch_files/SBATCH_get_unique_kmer_per_window.sh
6) Merge kmers present in sequences from a node and are specific on the node. Done using util.py -> get_unique_kmers_per_node(). This creates the nodekmers*txt files which lists all identified kmer markers for each of the node. These kmers are then used to run EM and threading algorithms.


