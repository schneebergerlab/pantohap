#!/bin/bash
#SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=15000mb
#SBATCH --time=12:00:00
#SBATCH -J QUickScript

module load samtools
chars=({A..J})
for i in 1 2 3 4; do
    if [ ! -f dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr${chrid}_hap${i}syri.out ]; then
  	  srun /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam \
  	  -n ${SLURM_CPUS_PER_TASK} -p dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr${chrid}_hap${i} \
  	  DM_chr${chrid}.fa ${chars[${SLURM_ARRAY_TASK_ID}]}_chr${chrid}_hap${i}.fa
  	else
  	  echo ; echo found dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr${chrid}_hap${i}syri.out ; echo
    fi
done

#module load samtools
##chars=({O..O})
##for i in 5 6 7 8; do
#srun /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools \
#  runsyri -alignment bam \
#	-n 10 -p dm_O_chr02_hap${SLURM_ARRAY_TASK_ID} \
#	DM_chr02.fa O_chr02_hap${SLURM_ARRAY_TASK_ID}.fa
##srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=1000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam -n 10 -p chr02_dm_A_hap${i} DM_chr02.fa chr02_hap${i}.fa
##done
##wait
#
#
