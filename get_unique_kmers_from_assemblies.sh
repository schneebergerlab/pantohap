#!/bin/bash
#SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=12:00:00
#SBATCH -J get_unique_kmers

indir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/assemblies/'
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
chars=({A..J})

for i in 1 2 3 4; do
  cd $cwd
  mkdir ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}; cd ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
  echo $(pwd)
  # Call kmers
  srun $meryl k=51 count threads=4 memory=20 output ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i} ${indir}/${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}_genome.fasta
  # Select unique kmers
  srun $meryl print equal-to 1 ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i} \
   | gzip > ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}.unique.k51.txt.gz
   # Tar and delete the database
  tar -cf ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}.tar ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
  rm -r ${chars[${SLURM_ARRAY_TASK_ID}]}_hap${i}
done