#!/bin/bash
#SBATCH --array=0
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=12:00:00
#SBATCH -J get_unique_kmers

chars=({A..J})
#for i in 5 6 7 8; do
for i in 5 ; do
	srun --exclusive --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
	  ../../tool/get_unique_kmer_per_window.sh ${chars[${SLURM_ARRAY_TASK_ID}]} $i 10 5 # threads and mem also parsed
done
wait