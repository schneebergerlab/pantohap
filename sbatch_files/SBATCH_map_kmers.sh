#!/bin/bash
#SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=1-00:00:00
#SBATCH -J map_kmers


#indir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/assemblies/'
indir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/data/potato/assemblies_v3/
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/kmer_size_51/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
bowtie2=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/bowtie2

chars=({A..J})
QUAL='~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

for h in 1 2 3 4; do
  {
    c=${chars[${SLURM_ARRAY_TASK_ID}]}
    sample=${c}_hap${h}
    cd ${cwd}/${sample}/ || exit

    srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
    $meryl print ${sample}_good 2> /dev/null \
    | awk -v QUAL="$QUAL" '{print "@kmer_"NR" "$2"\n"$1"\n+\n"QUAL}' \
    | pigz  -p ${SLURM_CPUS_PER_TASK} > ${sample}_good_kmers.fq.gz

    # Map kmers
    srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
    $bowtie2 -x ${indir}/hap4by12_${c}/${sample}_genome -U ${sample}_good_kmers.fq.gz -S ${sample}_good_kmers.sam --end-to-end --mp 10000 -p ${SLURM_CPUS_PER_TASK}

    # Create sorted bam file
    srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -O BAM ${sample}_good_kmers.sam > ${sample}_good_kmers.sorted.bam

    # Index bam
    srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
    samtools index -@ ${SLURM_CPUS_PER_TASK} ${sample}_good_kmers.sorted.bam

    # Split BAM files for individual chromosomes
    chrs=$(samtools view -H  ${sample}_good_kmers.sorted.bam | grep '@SQ' | cut -f2 | cut -d':' -f2 | grep -v 'scaffold')
    srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
    | {
      for c in ${chrs[@]}; do
        samtools view -@ ${SLURM_CPUS_PER_TASK} -O BAM ${sample}_good_kmers.sorted.bam $c > ${sample}_${c}_good_kmers.sorted.bam
        samtools index -@ ${SLURM_CPUS_PER_TASK} ${sample}_${c}_good_kmers.sorted.bam
      done
      }
  } &
done
wait
