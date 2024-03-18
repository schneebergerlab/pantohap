#!/bin/bash
genome=$1
hap=$2
cpu=$3
mem=$4
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
goodkmerdir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/

# Get 51-mers in each fasta sequence file
<<comment
ls syn_fasta_${genome}_hap${hap}*fasta \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} $meryl k=51 count threads=1 memory=5 output {} {}.fasta
comment

# Untar (if tar)
tar -xf syn_fasta_${genome}_hap${hap}.tar

# Select kmers that are present only once in the fasta and are good kmer in the sample
<<comment
tar -tf syn_fasta_${genome}_hap${hap}.fasta.tar.gz \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} bash -c "$meryl k=51 threads=1 memory=${mem} intersect [ equal-to 1 {} ] ${goodkmerdir}/${genome}_hap$((hap-4))/${genome}_hap$((hap-4))_good output {}_good " -- {}
comment

# Print selected kmers
tar -tf syn_fasta_${genome}_hap${hap}.fasta.tar.gz \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} bash -c "$meryl k=51 threads=1 memory=${mem} print {}_good  > {}.k51.good.txt " -- {}

<<comment
meryl union syn_fasta_${genome}_hap${hap}_bin_*_good output ${genome}_hap${hap}_good
comment

# Cleanup
# Archive meryl databases
tar -cf syn_fasta_${genome}_hap${hap}.tar $(ls -d syn_fasta_${genome}_hap${hap}_bin_*/)
rm -r $(ls -d syn_fasta_${genome}_hap${hap}_bin_*/)

# Archive fasta files
<<comment
tar -zcf syn_fasta_${genome}_hap${hap}.fasta.tar.gz syn_fasta_${genome}_hap${hap}*fasta
comment

# Archive selected kmers
## Currently, not used because util.get_unique_kmers_per_node is set up to require random access to kmer files
<<comment
tar -zcf syn_fasta_${genome}_hap${hap}.k51.good.txt.tar.gz syn_fasta_${genome}_hap${hap}_bin_*.k51.good.txt
rm syn_fasta_${genome}_hap${hap}_bin_*.k51.good.txt
comment

echo "done"
