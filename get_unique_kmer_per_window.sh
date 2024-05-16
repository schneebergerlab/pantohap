#!/bin/bash
genome=$1
hap=$2
k=$3
cpu=$4
mem=$5

meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
goodkmerdir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
synfastadir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/synfastas

# Extract fasta files from the tars
tar -xf ${synfastadir}/syn_fasta_${genome}_hap${hap}.fasta.tar.gz

# Get k-mers in each fasta sequence file
ls syn_fasta_${genome}_hap${hap}*fasta \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} $meryl k=${k} count threads=1 memory=5 output {} {}.fasta

<<comment
# Untar (if tar)
tar -xf syn_fasta_${genome}_hap${hap}.tar
comment

# Select kmers that are present only once in the fasta and are good kmer in the sample
ls syn_fasta_${genome}_hap${hap}*fasta \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} bash -c "$meryl k=${k} threads=1 memory=${mem} intersect {} ${genome}_hap$((hap-4))_good output {}_good " -- {}


# Print selected kmers
ls syn_fasta_${genome}_hap${hap}*fasta \
| sed 's/\.fasta//g' \
| xargs -n1 -P ${cpu} -I {} bash -c "$meryl k=${k} threads=1 memory=${mem} print {}_good  > {}.k51.good.txt " -- {}

<<comment
meryl union syn_fasta_${genome}_hap${hap}_bin_*_good output ${genome}_hap${hap}_good
comment

# Cleanup
# Archive meryl databases
tar -cf syn_fasta_${genome}_hap${hap}.tar $(ls -d syn_fasta_${genome}_hap${hap}_bin_*/)
rm -r $(ls -d syn_fasta_${genome}_hap${hap}_bin_*/)

# Archive fasta files
tar -zcf syn_fasta_${genome}_hap${hap}.fasta.tar.gz syn_fasta_${genome}_hap${hap}*fasta
rm syn_fasta_${genome}_hap${hap}*fasta

# Archive selected kmers
## Currently, not used because util.get_unique_kmers_per_node is set up to require random access to kmer files
<<comment
tar -zcf syn_fasta_${genome}_hap${hap}.k51.good.txt.tar.gz syn_fasta_${genome}_hap${hap}_bin_*.k51.good.txt
rm syn_fasta_${genome}_hap${hap}_bin_*.k51.good.txt
comment

echo "done"
