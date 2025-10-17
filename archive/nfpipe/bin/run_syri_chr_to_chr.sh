#!/bin/bash

chr=$1
sample=$2
hap=$3
threads=$4
cwd=$5
bam=$6

refgenome=${cwd}/data/DM_chr${chr}.fa
qrygenome=${cwd}/data/${sample}_chr${chr}_hap${hap}.fa
#minimap2 -ax asm5 --eqx -t $threads $refgenome $qrygenome >
#| samtools sort -O BAM - > dm_${sample}_chr${chr}_hap${hap}.bam
syri -c $bam -r $refgenome -q $qrygenome -F B --nc $threads --prefix dm_${sample}_chr${chr}_hap${hap} -f

