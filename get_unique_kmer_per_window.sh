#!/bin/bash
genome=$1
hap=$2
cpu=$3
mem=$4
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl

# Get 51-mers in each fasta sequence file
ls syn_fasta_${genome}_hap${hap}*fasta \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} $meryl k=51 count threads=1 memory=5 output {} {}.fasta

# Select kmers that are present only once in the fasta
ls syn_fasta_${genome}_hap${hap}*fasta \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} bash -c "$meryl print k=51 equal-to 1 threads=1 memory=5 {} > {}.k51.unique.txt " -- {}

# Cleanup
tar -cf syn_fasta_${genome}_hap${hap}.tar $(ls -d syn_fasta_${genome}_hap${hap}_bin_*/)
rm -r $(ls -d syn_fasta_${genome}_hap${hap}_bin_*/)
