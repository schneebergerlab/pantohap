#!/bin/bash
genome=$1
hap=$2
cpu=$3
mem=$4
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl

ls syn_fasta_${genome}_hap${hap}*fasta \
| xargs -n1 -P ${cpu} -I {} $meryl k=51 count threads=1 memory=5 output {} {}

#
#$meryl k=21 count threads=40 memory=100 output ${s}_R1.meryl ${s}_ql-trimmed-pair1.fastq.gz &
#    $meryl k=21 count threads=40 memory=100 output ${s}_R2.meryl ${s}_ql-trimmed-pair2.fastq.gz &
#    "
#done
#
## Combine Meryl databases
#$meryl union-sum output rp.meryl A_R1.meryl A_R2.meryl B_R1.meryl B_R2.meryl C_R1.meryl C_R2.meryl D_R1.meryl D_R2.meryl