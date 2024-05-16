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

cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
chars=({A..J})

for i in 5 6 7 8; do
  {
  for k in 21 31 41 51; do
    cd $cwd; cd kmer_size_${k}; cd ${chars[${SLURM_ARRAY_TASK_ID}]}_hap$((hap-4))

    srun --exclusive --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      ../../tool/get_unique_kmer_per_window.sh ${chars[${SLURM_ARRAY_TASK_ID}]} $i $k ${SLURM_CPUS_PER_TASK} 5 # threads and mem also parsed
  done
  } &
done
wait
