# <editor-fold desc='Step1: Generate haplotype graph'>
# Update SNPs genotype to include the information about the deletion or SR markers
# SNP_type - Genotype
# REF_allele in syntenic            -           0
# ALT_allele in syntenic            -           1
# Deleted                           -           2
# In inversion                      -           3
# In TDs (considered as deletion)   -           2


# Generate haplotype genotype file required for generating haplotype graph
get_haplotype_vcf.sh

# Currently the graph is generated using binning of reference genomes (util.py)
util.hapnodesfromvcffile()

# TODO: Consider using msyd hap-graph for this problem
run_msyd.sh
# </editor-fold>

# <editor-fold desc='Kmer counting'>

# For each genome, select kmers present only once in the assembly
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/
meryl=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/meryl
cd $cwd
sbatch ../../tool/sbatch_files/SBATCH_get_unique_kmers_from_assemblies_part1.sh


# Merge multimers
# SH
# for k in 21 31 41 51; do
for k in 51; do
    cd ${cwd}/kmer_size_${k}
    $meryl union  threads=40 memory=20 output all_genome_multi *hap*/*multi
done

sbatch ../../tool/sbatch_files/SBATCH_get_unique_kmers_from_assemblies_part2.sh

# Merge unimers
# SH
# for k in 21 31 41 51; do
for k in 51; do
    cd ${cwd}/kmer_size_${k}
    $meryl union output all_genome_good *hap*/*good &
done


# for each 100kb window, get syntenic sequence from query (40 hap) genomes.
cwd=/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/node_kmers/
util.get_node_query_sequence()

# Get kmers for each of the fasta file created using get_node_query_sequence()
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/node_kmers/
for c in {01..12}; do
    cd $cwd ; cd chr$c
    chrid=$c
    export chrid
    sbatch ../../../../tool/sbatch_files/SBATCH_get_unique_kmer_per_window.sh
done

# Merge Kmers from each collapsed node and then find marker kmers (that are present in that node only) for each node
#from functools import partial
#from multiprocessing import Pool
#cwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/kmer_analysis/'
#cwd = '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/'
#with Pool(processes=4) as pool:
#    pool.map(partial(get_unique_kmers_per_node, cwd=cwd), [21, 31, 41, 51])

# Haplotype graph filename
k=51
chrs = ["chr{:02d}".format(i) for i in range(1, 13)]
#for d in [0.1, 0.05, 0.01]:
for d in [0.1]:
#    hgf = 'haplotype_graph_{c}' + '_div{d}_2024_08_02.txt'.format(d=d)
    hgf = 'haplotype_graph_notd_{c}' + '_div{d}_2024_10_28.txt'.format(d=d)
    with Pool(processes=12) as pool:
        pool.map(partial(get_unique_kmers_per_node, hgf=hgf, k=k), chrs)


# Run EM to get ploidy levels for each node

# Get threads from the graph

# Plot threads
cwd=/home/ra98jam/d16/projects/potato_hap_example/results/threading/
indir=/home/ra98jam/d16/projects/potato_hap_example/data/

cd $cwd ; mkdir -p thread_plots
cd thread_plots
for c in {01..12}; do
    for g in WhiteRose Kenva RussetBurbank ; do
        {
        for m in bg th; do
#            ../threads_for_Manish_2024_08_28/${g}_chr${c}_div0.1.threads.tsv \
            python /home/ra98jam/d16/projects/potato_hap_example/tool/plot_hap_graph_2.py \
                ${indir}/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
                ../threading_to_report_final/threading_including_pericentromeres_final/${g}_chr${c}_div0.1.threads.tsv \
                ${g}_chr${c}_div0.1.threads.${m}.pdf -W 12 -H 4 --spread 5 --mode ${m}

            # python /home/ra98jam/d16/projects/potato_hap_example/tool/plot_hap_graph_2.py \
            #     ${indir}/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
            #     ../threads_plusCentromere_forManish_2024_08_27/${g}_chr${c}_div0.1.threads.tsv \
            #     ${g}_chr${c}_div0.1.threads.with_centro.${m}.pdf -W 12 -H 4 --spread 5 --mode ${m}
        done
        } &
    done
done

# Get sample-specific zoomed in view
c=06
g=WhiteRose
python /home/ra98jam/d16/projects/potato_hap_example/tool/plot_hap_graph_2.py \
    ${indir}/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
    ../threading_to_report_final/threading_including_pericentromeres_final/${g}_chr${c}_div0.1.threads.tsv \
    ${g}_chr${c}_div0.1.threads.haps.45_47.5mb.pdf -W 12 -H 4 --spread 5 --mode haps -s 45000000 -e 47500000 \
    --haplist A_hap1 --haplist A_hap2 --haplist A_hap3 --haplist A_hap4

python /home/ra98jam/d16/projects/potato_hap_example/tool/plot_hap_graph_2.py \
    ${indir}/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
    ../threading_to_report_final/threading_including_pericentromeres_final/${g}_chr${c}_div0.1.threads.tsv \
    ${g}_chr${c}_div0.1.threads.45_47.5mb.pdf -W 12 -H 4 --spread 5 --mode th -s 45000000 -e 47500000

python /home/ra98jam/d16/projects/potato_hap_example/tool/plot_hap_graph_2.py \
    ${indir}/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
    ../threading_to_report_final/threading_including_pericentromeres_final/${g}_chr${c}_div0.1.threads.tsv \
    ${g}_chr${c}_div0.1.threads.haps.45_50mb.pdf -W 12 -H 4 --spread 5 --mode haps -s 45000000 -e 50000000 \
    --haplist A_hap1 --haplist A_hap2 --haplist A_hap3 --haplist A_hap4

python /home/ra98jam/d16/projects/potato_hap_example/tool/plot_hap_graph_2.py \
    ${indir}/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
    ../threading_to_report_final/threading_including_pericentromeres_final/${g}_chr${c}_div0.1.threads.tsv \
    ${g}_chr${c}_div0.1.threads.45_50mb.pdf -W 12 -H 4 --spread 5 --mode th -s 45000000 -e 50000000


# Get pseudo assemblies
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/threading/thread_fastas/
cd $cwd
for c in {01..12}; do
    for g in WhiteRose Kenva RussetBurbank ; do
      {
#        python ../../../tool/get_fasta_seq_for_thread.py \
#            ../../../data/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
#            ../threads_forManish_26_08_24/${g}_chr${c}_div0.1.threads.tsv \
#            ../../kmer_analysis/node_kmers/fasta_len_in_nodes.csv \
#            ${g}_chr${c}_div0.1.no_centro.fa chr${c}
        # Select threads longer than 1 node
        awk '{if($2>1){print $0}}' ../threading_to_report_final/threading_including_pericentromeres_final/${g}_chr${c}_div0.1.threads.tsv > ${g}_chr${c}_div0.1.threads.t2.tsv
        # Get Contigs
        python ../../../tool/get_fasta_seq_for_thread.py \
            ../../../data/chr${c}/haplotype_graph_chr${c}_div0.1_2024_08_02.txt \
            ${g}_chr${c}_div0.1.threads.t2.tsv \
            ../../kmer_analysis/node_kmers/fasta_len_in_nodes.csv \
            ${g}_chr${c}_div0.1.with_centro.fa chr${c}
        } &
    done
done
grep -c '>' *fa | tr ':' ' ' | sed 's/_/\t/g' | cut -f1,2,5 > number_of_pseudo_contigs_per_chromosome.txt


# python /home/ra98jam/d16/projects/potato_hap_example/tool/get_fasta_seq_for_thread.py ../../../data/chr10/haplotype_graph_chr10_div0.1_2024_08_02.txt ../threads_forManish_26_08_24/WhiteRose_chr10_div0.1.threads.tsv ../../kmer_analysis/node_kmers/fasta_len_in_nodes.csv WhiteRose_chr10_div0.1.no_centr0.fa chr10


# Map and compare Russet Burbank
cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/threading/thread_fastas/
# merge RB contigs
cat RussetBurbank_chr*fa > RussetBurbank.pseudo_contigs.fa
# Merge RB haplotypes
indir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/assemblies/
cp ${indir}/R_hap*_genome.fasta .
cat R_hap*_genome.fasta > RussetBurbank.genome.fa

minimap2 -cx asm5 --eqx -t 10 RussetBurbank.genome.fa RussetBurbank.pseudo_contigs.fa > rb.contigs_to_genome.paf &
minimap2 -cx asm5 --eqx -t 10 RussetBurbank.pseudo_contigs.fa RussetBurbank.genome.fa > rb.genome_to_contigs.paf &


# Generate summary statistics
summary_plots_kmers_per_node()


# </editor-fold>
