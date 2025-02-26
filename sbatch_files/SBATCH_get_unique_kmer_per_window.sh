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

cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/node_kmers/
chars=({A..J})

#c=${chars[$((SLURM_ARRAY_TASK_ID/4))]}
#i=${is[$((SLURM_ARRAY_TASK_ID%4))]}

for i in 1 2 3 4; do
  {
#  for k in 21 31 41 51; do
  for k in 51; do
    cd $cwd; cd chr${chrid}; cd ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
    mkdir -p kmer_size_notd_${k} ; cd kmer_size_notd_${k}

    srun --exclusive --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
	    /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/tool/get_unique_kmer_per_window.sh \
      ${chars[${SLURM_ARRAY_TASK_ID}]} \
      $chrid \
      $i \
      $k \
      ${SLURM_CPUS_PER_TASK} \
      5 # threads and mem also parsed
  done
  } &
done
wait
