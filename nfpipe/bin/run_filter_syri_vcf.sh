#!/bin/bash

sed -i 's/v4.3/v4.2/' dm_${g}_chr${c}_hap${i}syri.vcf
grep -P -v 'TRANS|INV|DUP' dm_${g}_chr${c}_hap${i}syri.vcf > dm_${g}_chr${c}_hap${i}syri.nosr.vcf

finname=dm_${g}_chr${c}_hap${i}syri.nosr.snps
vcftools --vcf dm_${g}_chr${c}_hap${i}syri.nosr.vcf \
  --remove-indels --recode --recode-INFO-all \
  --chr chr${c} \
  --out $finname
sed -i "s/sample/${g}_hap${i}/" ${finname}.recode.vcf
bgzip ${finname}.recode.vcf
tabix ${finname}.recode.vcf.gz
bcftools sort -T . ${finname}.recode.vcf.gz -Oz -o ${finname}.sorted.vcf.gz
bcftools index  ${finname}.sorted.vcf.gz
