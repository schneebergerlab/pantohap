cwd='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/generate_hapgraph/'
datadir='/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/'

cd $cwd
# Setup genomes.csv file
rm genomes.csv
for c in {A..J}; do
  for i in 5 6 7 8; do
    echo -e ${c}_chr02_hap${i}' \t 'alignments/dm_${c}_chr02_hap${i}.bam' \t 'syriout/dm_${c}_chr02_hap${i}syri.out' \t 'syriout/dm_${c}_chr02_hap${i}syri.vcf' \t '${datadir}/${c}_chr02_hap${i}.fa >> genomes.csv
  done
done

# Run msyd
cd $cwd
module load samtools
sbatch --get-user-env \
  --clusters=biohpc_gen \
  --partition=biohpc_gen_normal \
  --ntasks=1 \
  --cpus-per-task=10 \
  --mem=20000mb \
  --time=12:00:00 \
  -J msyd_chr2 \
  /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/msyd call -i genomes.csv -c 10 -o chr02_all_hap.pff --realign -r ${datadir}/DM_chr02.fa

