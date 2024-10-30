"""
Use unique kmers from each genome and generate a haplotype-graph where
nodes are kmers and edges link kmers next to each other.
Kmers are sorted based on the topological order in the respective assemblies.
"""


# Get Kmers and map them to genomes
SBATCH_map_kmers.sh

# Assert that each kmer mapped perfectly
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/kmer_size_51/
for c in {A..J}; do
{
    for h in {1..4}; do
        sample=${c}_hap${h}
        cd ${cwd}/${sample}/
        n_kmers=$(meryl statistics ${sample}_good 2> /dev/null | grep unique | awk '{print $2}')
        n_mapped=$(samtools view ${sample}_good_kmers.sorted.bam | grep -c 'MD:Z:51')
        if [[ $n_kmers == $n_mapped ]]; then
            echo equal for $sample
        else
            echo not equal for $sample
        fi
    done
} &
done


# Get list of all good kmers for each chr
util.getkmerperchr()






cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/kmer_size_51/
cd $cwd
meryl threads=20 print union *hap*/*good 2> /dev/null | cut -f1 | pigz -p 20 > all_kmer.txt.gz


