#Kmer Pipeline for Matt

#1. Decompose your genomes into unique kmers and find their positions in their home genomes

#1a. Count 21mers in your reference genomes - you might want to make this higher (ie 31 should be overkill) so that you get more unique kmers later, but this might also expose some poor memory management in my scripts XD
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
nthread=10
mer=21
hash=80M
cd ${wd}/a1_kmer_decomposition/a1_genomes
for genome in chr06_hap21 chr06_hap22 chr06_hap23 chr06_hap24; do
    cd ${wd}/a1_kmer_decomposition/a1_genomes
    bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=32000]" -M 32000 -o jellyfish_${genome}.log -e jellyfish_${genome}.err "jellyfish count -m ${mer} -s ${hash} -t ${nthread} -C -o ${genome}_mer_counts.jf ../../d1_OtavaChr6s/${genome}.fasta"
    cd $wd
done
#produces a _mer_counts.jf file

# 1b. Now turn those counts into i) A histogram of kmer coverage (_histo.txt) and
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
for genome in chr06_hap21 chr06_hap22 chr06_hap23 chr06_hap24; do #
    cd ${wd}/a1_kmer_decomposition/a1_genomes
    jellyfish histo ${genome}_mer_counts.jf > ${genome}_mer_counts_histo.txt
    jellyfish dump ${genome}_mer_counts.jf > ${genome}_mer_counts_dumps.fa
    cd $wd
done

#example _histo.txt
#	1 40684573
#	2 2823782
#	3 821698
#	4 364780
#	5 202506
#	6 126809
#	7 86013
#	8 59805
#	9 45741
#	10 36984
#column 1 is the kmer frequency, column 2 is the number of kmers with that frequency in this genome

#example _counts_dumps.fa
#	>445
#	AAAAAAAAAAAAAAAAAAAAA
#	>1
#	TTATGCATTTTCCCATGAAAA
#	>1
#	GAAATTCACCAGTTGATCTGC
#	>1
#	GTCTTCAGAATTTGAAAGATA
#	>1
#	CCAGATATGGCTATACTTAAA
#Fasta header is the count for this kmer sequence, fasta sequence is the kmer sequence

#Now we want to find a set of kmers which are uniquely occuring in all genomes (For you maybe you only care about those which are unique across BOTH/ALL genomes, which we could select here, but I instead deal with later..)

#1c. Identify presence/absence of kmers and unique/repeated status
#Sort kmers from each genome - the .fa format above is a bit awkward, I wrote a script to invert the values and labels, then perform sorting, then invert it back to how it was.
scriptDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/r1_scripts
fastaFolder=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes
for genome in chr06_hap21 chr06_hap22 chr06_hap23 chr06_hap24; do
    #echo $genome
    bsub -q normal -R "span[hosts=1] rusage[mem=6000]" -M 6000 -o sort_${genome}.log -e sort_${genome}.err "bash ${scriptDir}/sort_jellyfishCounts_awk.sh ${fastaFolder}/${genome}_mer_counts_dumps.fa"
    #using the awk version of this script going forward, because the bash-only version was toooo slowwww (5 hrs)
done

#1d. Iterate through sorted kmers to find those that are unique and not multiple-occurring
#Identify set of kmers which are either unique or missing in all genomes
#output both a table with presence/absence of all kmers, and a fasta file of those kmers to be used for counting in short reads/single-cell later.
fastaFolder=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes
outDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r1_merged_kmer_counts_forUnique
python ${scriptDir}/findUniqueKmers_acrossGenomes.py ${fastaFolder}/chr06_hap21_mer_counts_dumps.sorted.fa,${fastaFolder}/chr06_hap22_mer_counts_dumps.sorted.fa

#1e.  identify locations of each kmer in the respective genomes
#Map kmers back to genomes with bwa
#i) - index Haplotypes (already done)
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/d1_OtavaChr6s
cd ${wd}
for genome in chr06_hap21 chr06_hap22 chr06_hap23 chr06_hap24; do
bsub -m "hpc001 hpc002 hpc003 hpc004 hpc005 hpc006" -q normal -R "rusage[mem=5000]" -R "span[hosts=1]" -M 5000 -o index.log -e index.err "bwa index ${wd}/${genome}.fasta"
done

#ii) - Align with BWA MEM
#This gives us a sam file with mapping of reads in the same order as merged kmers file
#You might want BAMs/CRAMS, I just used SAMs for testing so that I could sanity check each step. ANd my later tools expect to read in SAMs - early days on all this.
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
cd ${wd}
for hap in chr06_hap21 chr06_hap22 chr06_hap23 chr06_hap24; do
    kmers=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r1_merged_kmer_counts_forUnique/Otava_Chr06.fa
    kmerSet=$(basename "${kmers}" .fa)
    genome=${wd}/d1_OtavaChr6s/${hap}.fasta
    outDir=${wd}/a1_kmer_decomposition/a1_genomes/r2_map_kmers_backToGenomes
    cd $outDir
    bsub -o ${hap}_bwa_mem.log -e ${hap}_bwa_mem.err -q multicore20 -n 8 -R "rusage[mem=16000]" -R "span[hosts=1]" -M 16000 "bwa mem -k 21 -T 21 -a -c 1 -t 8 -M -Y ${genome} ${kmers} > ${outDir}/${kmerSet}_kmers_to_${hap}.bwaAligned.sam"
    #| samtools view -h -f 0 -F 256 -o ${outDir}/${kmerSet}_kmers_to_${hap}.bwaAligned.sam - "
    cd ${wd}
done

#Example SAM output at this stage
#
#	@SQ	SN:chr06_hap22	LN:57327765
#	@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:/opt/share/software/packages/bwa-0.7.17/bin/bwa mem -k 21 -T 21 -a -c 1 -t 8 -M -Y /netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/d1_OtavaChr6s/chr06_hap22.fasta /netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r1_merged_kmer_counts_forUnique/Otava_Chr06.fa
#	1	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAAACGG	*	AS:i:0	XS:i:0
#	2	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAAATGT	*	AS:i:0	XS:i:0
#	3	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAACAGC	*	AS:i:0	XS:i:0
#	4	16	chr06_hap22	47290450	0	21M	*	0	0	GATGTTTTTTTTTTTTTTTTT	*	NM:i:0MD:Z:21	AS:i:21	XS:i:0
#	5	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAACCGA	*	AS:i:0	XS:i:0
#	6	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAACCGG	*	AS:i:0	XS:i:0
#	7	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAACGAC	*	AS:i:0	XS:i:0
#	8	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAACGCC	*	AS:i:0	XS:i:0
#	9	4	*	0	0	*	*	0	0	AAAAAAAAAAAAAAAAACGGA	*	AS:i:0	XS:i:0
#	10	16	chr06_hap22	48098729	0	21M	*	0	0	ACCGTTTTTTTTTTTTTTTTT	*	NM:i:0MD:Z:21	AS:i:21	XS:i:0


#Okay, Now we know what kmers we care about to use as markers, and we know where they occur.

#2a
#This is basically the same procedure as 1a-b, but there is a lot more to count in a fasta file. Jellyfish can get around this with bloom filters, but the counts can be innacurate around low read coverages (which is most of what you have)
#But we have a secret weapon. We already know which kmers we want to count! (from 1c) So we can do it pretty fast, here I dind't even use the cluster.
#This step I guess you woulf want to repeat
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
readset=chr06_hap21_20x
cd ${wd}/a1_kmer_decomposition/a2_shortReads
mkdir $readset
cd $readset
merged_fasta=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r1_merged_kmer_counts_forUnique/Otava_Chr06.fa
r1=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/d2_simulated_WGS/chr06_hap21_20x/chr06_hap21.20x.r1.fq.gz
r2=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/d2_simulated_WGS/chr06_hap21_20x/chr06_hap21.20x.r2.fq.gz
jellyfish count -m 21 -s 240M -C -o ${readset}_21mer_counts.mergedOnly.jf --if ${merged_fasta} <(zcat ${r1}) <(zcat ${r2})
#--lower-count 2 #Try without this, so that we still get 0 counts?
#all runs pretty quickly
#   #2a   #Generate kmer count histogram
jellyfish histo ${readset}_21mer_counts.mergedOnly.jf> ${readset}_21mer_counts.mergedOnly.histo.txt
jellyfish dump ${readset}_21mer_counts.mergedOnly.jf > ${readset}_21mer_counts.mergedOnly.dumps.fa

#sort the output fasta
bsub -q normal -R "span[hosts=1] rusage[mem=6000]" -M 6000 -o sort_${readset}.log -e sort_${readset}.err "bash ${scriptDir}/sort_jellyfishCounts_awk.sh  ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.dumps.fa"

####
#Here I would then model the kmer count distribution. But I guess that you only have counts of 0,1,2 so can skip this
#wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
#scriptDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/r1_scripts
#readset=chr06_hap21_20x
#python ${scriptDir}/fit_peaks_v4.py fit -H ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.histo.txt -p 1 -e 20

#Then I would store the count and presence/absence information into the SAM file, as custom tags: so that my positional, count, and presence/absence information is all together.
#But further below I've provided a modified version which doesn't need this probabilty model stuff above
#wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
#scriptDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/r1_scripts
#readset=chr06_hap21_20x
#python ${scriptDir}/fit_peaks_v4.py prob -M ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.histo.txt.gaussianFit.tsv -C ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.dumps.sorted.fa
#python ${scriptDir}/fit_peaks_v4.py prob_display -M ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.histo.txt.gaussianFit.tsv -x 0 -X 200
#samDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r2_map_kmers_backToGenomes
#SAMList=${samDir}/Otava_Chr06_kmers_to_chr06_hap21.bwaAligned.sam,${samDir}/Otava_Chr06_kmers_to_chr06_hap22.bwaAligned.sam,${samDir}/Otava_Chr06_kmers_to_chr06_hap23.bwaAligned.sam,${samDir}/Otava_Chr06_kmers_to_chr06_hap24.bwaAligned.sam
#ProbFile=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a2_shortReads/chr06_hap21_20x/chr06_hap21_20x_21mer_counts.mergedOnly.dumps.sorted.fa.gaussianFit.ploidy_1.tsv
#PresenceFile=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r1_merged_kmer_counts_forUnique/Otava_Chr06.PresenceMatrix.tsv
#
#python ${scriptDir}/fit_peaks_v4.py write_sam_tags -P ${ProbFile} -S ${SAMList} -o ${samDir}/taggedSAMs -a ${PresenceFile}

###Not working for your purpose
#Assumes 1-to-1 match of the lines of the SAM file made in step 1d and the
wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica
scriptDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/r1_scripts
readset=chr06_hap21_20x
python ${scriptDir}/fit_peaks_v4.py prob -M ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.histo.txt.gaussianFit.tsv -C ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.dumps.sorted.fa
python ${scriptDir}/fit_peaks_v4.py prob_display -M ${wd}/a1_kmer_decomposition/a2_shortReads/${readset}/${readset}_21mer_counts.mergedOnly.histo.txt.gaussianFit.tsv -x 0 -X 200
samDir=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r2_map_kmers_backToGenomes
SAMList=${samDir}/Otava_Chr06_kmers_to_chr06_hap21.bwaAligned.sam,${samDir}/Otava_Chr06_kmers_to_chr06_hap22.bwaAligned.sam,${samDir}/Otava_Chr06_kmers_to_chr06_hap23.bwaAligned.sam,${samDir}/Otava_Chr06_kmers_to_chr06_hap24.bwaAligned.sam
ProbFile=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a2_shortReads/chr06_hap21_20x/chr06_hap21_20x_21mer_counts.mergedOnly.dumps.sorted.fa.gaussianFit.ploidy_1.tsv
PresenceFile=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_AncestralReferenceAlignment/a5_KaptainMERica/a1_kmer_decomposition/a1_genomes/r1_merged_kmer_counts_forUnique/Otava_Chr06.PresenceMatrix.tsv

python ${scriptDir}/WriteSamTags.py write_sam_tags -P ${ProbFile} -S ${SAMList} -o ${samDir}/taggedSAMs -a ${PresenceFile}