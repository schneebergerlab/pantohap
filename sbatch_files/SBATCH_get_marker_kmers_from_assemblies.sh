#!/bin/bash
#SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=12:00:00
#SBATCH -J get_unique_kmers

#indir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/assemblies/'
indir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/data/potato/assemblies_v3/
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
chars=({A..J})
kmer=11
for i in 1 2 3 4; do
  {
#    for k in 21 31 41 51; do
    for k in ${kmer}; do
      # Create folder corresponding to specific K-mer for each haplotype
      cd $cwd
      mkdir -p marker_kmers_${k}; cd marker_kmers_${k}
      mkdir -p ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}; cd ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
      echo $(pwd)

      # Count kmers
      srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      $meryl k=${k} count threads=4 memory=20 output ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i} ${indir}/hap4by12_${chars[${SLURM_ARRAY_TASK_ID}]}/${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}_genome.fasta

      # Get multimeric kmers
      srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      $meryl greater-than 1 output ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}_multi ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}

    done
  } &
done
wait
#
## Merge all Kmers
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/marker_kmers_${kmer}/
cd $cwd
ls -d  *_hap*/*_hap* | grep -v 'multi' | xargs meryl union threads=20 memory=40 output all_genome_kmers

# Get unique kmers
meryl equal-to 1 threads=20 memory=40 output all_genome_kmers_uni all_genome_kmers

# Get Genome specific unique kmers
for s in {A..J}; do
  echo $s
  for i in 1 2 3 4; do
    cd $cwd
    cd ${s}_hap${i}
    meryl threads=1 memory=2 output ${s}_hap${i}_uni intersect [equal-to 1 ${s}_hap${i}] ../all_genome_kmers_uni &
  done
done

# Align kmers to genome assemblies
indir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/data/potato/assemblies_v3/
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/marker_kmers_${kmer}/

for s in {A..J}; do
  echo $s
  for i in 1 2 3 4; do
    {
    cd $cwd
    cd ${s}_hap${i}
    # Kmers to fasta
#    meryl print ${s}_hap${i}_uni 2> /dev/null | awk '{print "@"NR"\n"$1"\n+\n~~~~~~~~~~~~~~~~~~~~~"}' | bgzip > ${s}_hap${i}_uni.fq.gz
#    bowtie2-build ${indir}/hap4by12_${s}/${s}_hap${i}_genome.fasta ${indir}/hap4by12_${s}/${s}_hap${i}_genome
    bowtie2 --end-to-end \
      --very-sensitive \
      --threads 1 \
      -x ${indir}/hap4by12_${s}/${s}_hap${i}_genome  \
      -U ${s}_hap${i}_uni.fq.gz \
      --rg-id ${s}_hap${i}_uni \
    | samtools sort -@1 -O BAM - \
    > ${s}_hap${i}_uni.sorted.bam
    samtools index -@1 ${s}_hap${i}_uni.sorted.bam
    } &
  done
done