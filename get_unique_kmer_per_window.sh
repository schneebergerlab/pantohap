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
for i in 5 6 7 8; do
	srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=1000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl -h

#	/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam \
#	-n 10 -p dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i} \
#	DM_chr02.fa ${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i}.fa
#srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=1000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam -n 10 -p chr02_dm_A_hap${i} DM_chr02.fa chr02_hap${i}.fa
done
wait
#wait





#
#$meryl k=21 count threads=40 memory=100 output ${s}_R1.meryl ${s}_ql-trimmed-pair1.fastq.gz &
#    $meryl k=21 count threads=40 memory=100 output ${s}_R2.meryl ${s}_ql-trimmed-pair2.fastq.gz &
#    "
#done
#
## Combine Meryl databases
#$meryl union-sum output rp.meryl A_R1.meryl A_R2.meryl B_R1.meryl B_R2.meryl C_R1.meryl C_R2.meryl D_R1.meryl D_R2.meryl
#
