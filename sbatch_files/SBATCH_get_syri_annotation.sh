#!/bin/bash
#SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=12:00:00
#SBATCH -J get_syri_anno

chars=({A..J})

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr${chrid}

function get_anno () {
  tabix ${1}.out.bed.gz -R $fname \
  | sort -k1,1 -k2,2n \
  | uniq > ${1}.out.bed.snp_anno.txt
}

export -f get_anno

for i in 1 2 3 4; do
  {
    get_anno dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr${chrid}_hap${i}syri

  } &
done
wait
