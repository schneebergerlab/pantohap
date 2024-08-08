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

for i in 1 2 3 4; do
  {
#    for k in 21 31 41 51; do
    for k in 51; do
      # Create folder corresponding to specific K-mer for each haplotype
      cd $cwd
      mkdir -p kmer_size_${k}; cd kmer_size_${k}
      mkdir -p ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}; cd ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
      echo $(pwd)

      # Count kmers
      srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      $meryl k=${k} count threads=4 memory=20 output ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i} ${indir}/hap4by12_${chars[${SLURM_ARRAY_TASK_ID}]}/${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}_genome.fasta

      # Get unique kmers
      srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      $meryl equal-to 1 output ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}_uni ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}

      # Get multimeric kmers
      srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      $meryl greater-than 1 output ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}_multi ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}

      # Select unique kmers
    #  srun $meryl print equal-to 1 ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i} \
    #   | gzip > ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}.unique.k51.txt.gz
       # Tar and delete the database
      tar -cf ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}.tar ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
      rm -r ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
    done
  } &
done

wait
