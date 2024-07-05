#############################################################
# Step 1: Get chromosome fasta files
#############################################################
#A haplotype graph would represent one chromosome, so for each
# chromosome we create a different sets of input files.
#Lets say we are processing chr1, then we would need one fasta
#file consisting of DM:chr1 (lets say named: DM_chr01.fa) and N
#other fasta files (currently 40 for 40 haplotypes that we have)
# consisting of chromosome each lets say names (A_chr01_hap5.fa,
# A_chr01_hap6.fa etc).

#############################################################
# Step 2: Run syri
#############################################################

# Normal syri run between each DM and each haplotype

# My slurm script for chr02:
chars=({A..J})
for i in 5 6 7 8; do
	srun /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam \
	-n 10 -p dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i} \
	DM_chr02.fa ${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i}.fa
done

# ChatGPT translated equivalent script for LSF

#!/bin/bash
# Load the required modules or activate the environment if needed
# module load minimap2 samtools
# Define the array of characters
chars=({A..J})
for i in 5 6 7 8; do
  bsub -n 10 -R "span[hosts=1]" "hometools runsyri -alignment bam -n 10 -p dm_${chars[${LSB_JOBINDEX}]}_chr02_hap${i} DM_chr02.fa ${chars[${LSB_JOBINDEX}]}_chr02_hap${i}.fa"
done

#############################################################
# Step 3: Index syri.out files
#############################################################
# This is required for doing random access of syri.out files
# Index syri out as well
for c in {A..J}; do
    echo $c
    for h in {5..8}; do
        hometools syriidx dm_${c}_chr02_hap${h}syri.out &
    done
    wait
done


#############################################################
# Step 4: Filter syri.vcf files to get syntenic SNPs in all samples
#############################################################

# Filter out SRs and all variants within SR regions from the VCF file for each haplotype
for c in {A..J}; do
    echo $c
    for i in {5..8}; do
        if [ -f dm_${c}_chr02_hap${i}syri.vcf ]; then
            sed -i 's/v4.3/v4.2/' dm_${c}_chr02_hap${i}syri.vcf
            grep -P -v 'TRANS|INV|DUP' dm_${c}_chr02_hap${i}syri.vcf > dm_${c}_chr02_hap${i}syri.nosr.vcf
        fi
    done
done


# Select SNPs from the syri outputs and create sorted and indexed vcf files
for c in {A..J}; do
    echo $c
    for i in {5..8}; do
        # Check that the file exists
        if [ -f dm_${c}_chr02_hap${i}syri.vcf ]; then
            finname=dm_${c}_chr02_hap${i}syri.nosr.snps
            vcftools --vcf dm_${c}_chr02_hap${i}syri.nosr.vcf \
                --remove-indels --recode --recode-INFO-all \
                --chr chr02 \
                --out $finname
            sed -i "s/sample/${c}_hap${i}/" ${finname}.recode.vcf
            bgzip ${finname}.recode.vcf
            tabix ${finname}.recode.vcf.gz
            bcftools sort -T . ${finname}.recode.vcf.gz -Oz -o ${finname}.sorted.vcf.gz
            bcftools index  ${finname}.sorted.vcf.gz
        fi
    done
done

# Merge the 40 VCF files
bcftools merge dm*syri.nosr.snps.sorted.vcf.gz -Oz -o dm_all_sample_chr2.syri.nosr.snps.merged.vcf.gz
bcftools index dm_all_sample_chr2.syri.nosr.snps.merged.vcf.gz


#############################################################
# Step 5: Get genotype table
#############################################################

# Get a genotype table for the selected SNPs (Contains 1 for alt allele and '.' otherwise)
zgrep -v '^##' dm_all_sample_chr2.syri.nosr.snps.merged.vcf.gz \
| cut --complement -f 3,8 \
| grep -v ','
> dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt

# The following command does not seem to be necessary. I keep
# it here in case the ".regions" files are used somewhere.
tail -n +2 dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt | awk '{print $1"\t"$2"\t"$2}' > dm_all_sample_chr2.syri.nosr.snps.merged.vcf.regions


#############################################################
# Step 6: Get syri annotation at the selected SNP position
#############################################################

fname=dm_all_sample_chr2.syri.nosr.snps.merged.vcf.regions_for_tabix.txt
awk 'NR>1 { print $1"\t"$2-1"\t"$2+1}' dm_all_sample_chr2.syri.nosr.snps.merged.vcf.txt > $fname
for c in {A..J}; do
    for h in {5..8}; do
        tabix dm_${c}_chr02_hap${h}syri.out.bed.gz -R $fname | sort | uniq > dm_${c}_chr02_hap${h}syri.out.bed.snp_anno.txt &
    done
    wait
done

#############################################################
# Step 7: Change the genotype table to include deletion markers
#############################################################
# Check each sample at each SNP position and change the genotypes
# based on the annotation of the region in syri.

# SNP_type                          -        Genotype
# ---------------------------------------------------
# REF_allele in syntenic            -           0
# ALT_allele in syntenic            -           1
# In inversion                      -           2
# Deleted                           -           3
# In TDs (considered as deletion)   -           3

Run the python code inside util.updatesnpgenotypes (util.py)

# The output file dm_all_sample_chr2.syri.nosr.snps.merged.vcf.del_markers.txt
# is used as input to generate the graph.