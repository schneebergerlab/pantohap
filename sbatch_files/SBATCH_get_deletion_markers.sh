#!/bin/bash
###SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=12:00:00
#SBATCH -J get_syri_anno


cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr${chrid}
cd $cwd

srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
python ../../tool/get_deletion_markers.py $cwd chr${chrid} ${SLURM_CPUS_PER_TASK}

