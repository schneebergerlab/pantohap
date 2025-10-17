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
indir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/data/potato/assemblies_v3/
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/
for c in {01..12}; do
  cd $cwd; mkdir -p chr$c ; cd chr$c
  {
  for g in {A..J}; do
    for h in {1..4}; do
      chrid=$(( ((10#${c}-1)*4) + ${h} ))
      hometools getchr ${indir}/hap4by12_${g}/${g}_hap${h}_genome.fasta --chrs chr${c}_hap${chrid}_v3 -o ${g}_chr${c}_hap${h}.fa
#      echo $c $g $h
#      echo $c $g $h  $(( ((10#${c}-1)*4) + ${h} ))
    done
  done
  } &
done

# Get DM chromosomes
for c in {01..12}; do
  cd $cwd; cd chr$c
  hometools getchr ../../../../../pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz --chrs chr${c} -o DM_chr${c}.fa
#../../../../pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz
done

#############################################################
# Step 2: Run syri
#############################################################

# Normal syri run between DM and each haplotype

for c in {01..12}; do
  cd $cwd; cd chr$c
  chrid=$c
  export chrid
  sbatch ../../tool/sbatch_files/SBATCH_runsyri.sh
done

# These jobs (specifically minimap2 was randomly but consistently crashing on SLURM, so I did minimap2/syri runs locally on the local computer)
cwd=/home/ra98jam/d16/projects/potato_hap_example/data/
for c in {01..12}; do
  cd $cwd; cd chr$c
  {
  for g in {A..J}; do
    for i in 1 2 3 4; do
      if [ ! -f dm_${g}_chr${c}_hap${i}syri.summary ]; then
        echo $c $g $i
        hometools runsyri -alignment bam \
        -n 1 -p dm_${g}_chr${c}_hap${i} \
        DM_chr${c}.fa ${g}_chr${c}_hap${i}.fa 2> nohup.out
      else
        echo ; echo found dm_${g}_chr${c}_hap${i}syri.out ; echo
      fi
    done
  done
} &
done


#############################################################
# Step 3: Index syri.out files
#############################################################
# This is required for doing random access of syri.out files
# Index syri out as well
for c in {01..12}; do
  cd $cwd; cd chr$c
  for g in {A..J}; do
      for h in {1..4}; do
#        if [ ! -f dm_${g}_chr${c}_hap${h}syri.out.bed.gz ]; then
        hometools syriidx dm_${g}_chr${c}_hap${h}syri.out &
#        fi
      done
  done
  wait
done

# Index syri out for only SYN regions. This is used for fetching node_query_sequence
for c in {01..12}; do
  cd $cwd; cd chr$c
  for g in {A..J}; do
      for h in {1..4}; do
        {
#        if [ ! -f dm_${g}_chr${c}_hap${h}syri.out.bed.gz ]; then
        zgrep 'SYN' dm_${g}_chr${c}_hap${h}syri.out.bed.gz | bgzip > dm_${g}_chr${c}_hap${h}syri.out.syn.bed.gz
        tabix -fp bed dm_${g}_chr${c}_hap${h}syri.out.syn.bed.gz
#        fi
        } &
      done
  done
  wait
done

#############################################################
# Step 4: Filter syri.vcf files to get syntenic SNPs in all samples
#############################################################
# Syri runs done locally used older version of syri that does not output sample column in the VCF file. Adding it manually.
for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  for g in {A..J}; do
    for i in {1..4}; do
        if ! grep -q sample dm_${g}_chr${c}_hap${i}syri.vcf ; then
          sed -i '/^#CHROM/s/$/\tFORMAT\tsample/;/^chr/s/$/\tGT\t1/' dm_${g}_chr${c}_hap${i}syri.vcf &
        fi
    done
  done
wait
done


for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  for g in {A..J}; do
    for i in {1..4}; do
        if ! grep -q 'ID=GT' dm_${g}_chr${c}_hap${i}syri.vcf ; then
          sed -i '/^#CHROM/s/^/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n/' dm_${g}_chr${c}_hap${i}syri.vcf &
        fi
    done
  done
wait
done
#
#sed -i '/^#CHROM/s/^/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n/' BKP_dm_A_chr01_hap2syri.vcf
#sed '/^#CHROM/s/^/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n/' BKP_dm_A_chr01_hap2syri.vcf > tmp.txt
#


# Filter out SRs and all variants within SR regions from the VCF file for each haplotype
for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  for g in {A..J}; do
      for i in {1..4}; do
        {
          if [ -f dm_${g}_chr${c}_hap${i}syri.vcf ]; then
            sed -i 's/v4.3/v4.2/' dm_${g}_chr${c}_hap${i}syri.vcf
#            grep -P -v 'TRANS|INV|DUP' dm_${g}_chr${c}_hap${i}syri.vcf > dm_${g}_chr${c}_hap${i}syri.nosr.vcf
            grep -P -v 'TRANS|INVTR|DUP|INVDP' dm_${g}_chr${c}_hap${i}syri.vcf > dm_${g}_chr${c}_hap${i}syri.notd.vcf
          else
            echo missing dm_${g}_chr${c}_hap${i}syri.vcf
          fi
        } &
      done
  done
  wait
done


# Select SNPs from the syri outputs and create sorted and indexed vcf files
for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  for g in {A..J}; do
    for i in {1..4}; do
      # Check that the file exists
      {
      if [ -f dm_${g}_chr${c}_hap${i}syri.vcf ]; then
          echo "starting"
#          finname=dm_${g}_chr${c}_hap${i}syri.nosr.snps
#          vcftools --vcf dm_${g}_chr${c}_hap${i}syri.nosr.vcf \
          finname=dm_${g}_chr${c}_hap${i}syri.notd.snps
          vcftools --vcf dm_${g}_chr${c}_hap${i}syri.notd.vcf \
              --remove-indels --recode --recode-INFO-all \
              --chr chr${c} \
              --out $finname
          sed -i "s/sample/${g}_hap${i}/" ${finname}.recode.vcf
          bgzip ${finname}.recode.vcf
          tabix ${finname}.recode.vcf.gz
          bcftools sort -T . ${finname}.recode.vcf.gz -Oz -o ${finname}.sorted.vcf.gz
          bcftools index  ${finname}.sorted.vcf.gz
      fi
      } &
    done
  done
  wait
done

# Merge the 40 VCF files
for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  {
#  bcftools merge dm*syri.nosr.snps.sorted.vcf.gz -Oz -o dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.gz
#  bcftools index dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.gz
  bcftools merge dm*syri.notd.snps.sorted.vcf.gz -Oz -o dm_all_sample_chr${c}.syri.notd.snps.merged.vcf.gz
  bcftools index dm_all_sample_chr${c}.syri.notd.snps.merged.vcf.gz
  } &
done

#############################################################
# Step 5: Get genotype table
#############################################################

for c in {01..12}; do
  {
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  # Get a genotype table for the selected SNPs (Contains 1 for alt allele and '.' otherwise)
#  zgrep -v '^##' dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.gz \
#  | cut --complement -f 3,8 \
#  | grep -v ',' \
#  > dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.txt
  zgrep -v '^##' dm_all_sample_chr${c}.syri.notd.snps.merged.vcf.gz \
  | cut --complement -f 3,8 \
  | grep -v ',' \
  > dm_all_sample_chr${c}.syri.notd.snps.merged.vcf.txt

  # The following command does not seem to be necessary. I keep
  # it here in case the ".regions" files are used somewhere.
  tail -n +2 dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.txt | awk '{print $1"\t"$2"\t"$2}' > dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.regions
  } &
done

#############################################################
# Step 6: Get syri annotation at the selected SNP position
#############################################################
for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
#  fname=dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.regions_for_tabix.txt
#  awk 'NR>1 { print $1"\t"$2-1"\t"$2}' dm_all_sample_chr${c}.syri.nosr.snps.merged.vcf.txt > $fname
  fname=dm_all_sample_chr${c}.syri.notd.snps.merged.vcf.regions_for_tabix.txt
  awk 'NR>1 { print $1"\t"$2-1"\t"$2}' dm_all_sample_chr${c}.syri.notd.snps.merged.vcf.txt > $fname
  chrid=$c
  export chrid
  export fname
  sbatch ../../tool/sbatch_files/SBATCH_get_syri_annotation.sh
done


#    # With Otava
#    fname=dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.regions_for_tabix.txt
#    awk 'NR>1 { print $1"\t"$2-1"\t"$2+1}' dm_all_sample_chr2.with_Otava.syri.nosr.snps.merged.vcf.txt > $fname
#    for c in {A..J} O; do
#        {
#        for h in {5..8}; do
#            tabix dm_${c}_chr02_hap${h}syri.out.bed.gz -R $fname | sort | uniq > dm_${c}_chr02_hap${h}.with_Otava.syri.out.bed.snp_anno.txt
#        done
#        } &
#    done


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

for c in {01..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  chrid=$c
  export chrid
  sbatch ../../tool/sbatch_files/SBATCH_get_deletion_markers.sh
done


# The output file dm_all_sample_chr2.syri.nosr.snps.merged.vcf.del_markers.txt
# is used as input to generate the graph.