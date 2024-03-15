#!/bin/bash
###SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=12:00:00
#SBATCH -J get_unique_kmers

indir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/assemblies/'
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
chars=({A..J})

cd $cwd
echo $(pwd)
# Merge multimers
srun $meryl union output all_genome_multi */*multi

# Get unique kmers after removing bad kmers
for c in {A..J}; do
  for i in 1 2 3 4; do
    cd ${cwd}/${c}_hap${i}
    srun  --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=5000 $meryl difference  threads=1 memory=5 output ${c}_hap${i}_good ${c}_hap${i}_uni ../all_genome_multi &
  done
  wait
done

# Merge unimers
srun $meryl union output all_genome_good */*good