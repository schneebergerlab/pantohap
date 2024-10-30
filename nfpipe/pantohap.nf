//chrids = Channel.of('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12')
chrids = Channel.of('01', '02')
//samples = Channel.of("A", "B", "C", "D", "E", "F", "G","H","I", "J")
samples = Channel.of("A", "B")
//haps = Channel.of(1, 2 , 3, 4)
haps = Channel.of(1, 2)
indexes = chrids.combine(samples.combine(haps))
references = Channel.fromPath("data/DM*.fa")
chromosomes = Channel.fromPath("data/*chr*_hap*.fa")

process CREATECHROMINDICES {
	maxForks 20
	publishDir "data/"

	input:
	path fin
	
	output:
	path "${fin}.fai"
	
	script:
	"""
	samtools faidx $fin
	"""
}


process MINIMAP2_CHR_TO_CHR {
    fair true
    label 'MM2_ARRAY'

	input:
	tuple val(chr), val(sample), val(hap)

	output:
	path "dm_${sample}_chr${chr}_hap${hap}.sam"

	script:
	"""
	$projectDir/bin/run_minimap2_chr_to_chr.sh $chr $sample $hap ${task.cpus} $projectDir
	"""
}


process SAM_TO_BAM {
    fair true
    label 'S2B_ARRAY'

	input:
	tuple val(chr), val(sample), val(hap)
	path sam

	output:
	path "dm_${sample}_chr${chr}_hap${hap}.bam"

	script:
	"""
	samtools sort -O BAM -o dm_${sample}_chr${chr}_hap${hap}.bam $sam
	"""
}


process RUN_SYRI {
    fair true
    label 'SYRI_ARRAY'

    publishDir "results/syriout", pattern: "dm_${sample}_chr${chr}_hap${hap}{syri.vcf,syri.out}"

	input:
	tuple val(chr), val(sample), val(hap)
	path bam

	output:
	path "dm_${sample}_chr${chr}_hap${hap}syri.out"
	path "dm_${sample}_chr${chr}_hap${hap}syri.vcf"

	script:
	"""
	echo $chr $sample $hap $bam
	$projectDir/bin/run_syri_chr_to_chr.sh $chr $sample $hap ${task.cpus} $projectDir $bam
	"""
}


process SYRI_INDEX {
    fair true
    label 'SYRI_INDEX'

	input:
	path sout

	output:
	tuple path("${sout}.bed.gz"), path("${sout}.bed.gz.tbi")

	script:
	"""
    hometools syriidx ${sout}
    #zgrep 'SYN' ${sout}.bed.gz | bgzip > ${sout}.syn.bed.gz
    #tabix -fp bed ${sout}.syn.bed.gz
	"""

}


process SYRI_FILTER {
    fair true
    label 'syri_filter'

    input:
    tuple val(chr), val(sample), val(hap)
    path svcf

	output:
	tuple val(chr), path("${svcf}.notd.snps.sorted.vcf.gz"), path("${svcf}.notd.snps.sorted.vcf.gz.csi")

    script:
    """
    sed -i 's/v4.3/v4.2/' $svcf
    # grep -P -v 'TRANS|INV|DUP' ${svcf} > ${svcf}.nosr.vcf

    grep -P -v 'TRANS|INVTR|DUP|INVDP' ${svcf} > ${svcf}.notd.vcf


    #vcftools --vcf ${svcf}.nosr.vcf \
    #  --remove-indels --recode --recode-INFO-all \
    #  --chr chr${chr} \
    #  --out ${svcf}.nosr.snps
    #sed -i "s/sample/${sample}_hap${hap}/" ${svcf}.nosr.snps.recode.vcf
    #bgzip ${svcf}.nosr.snps.recode.vcf
    #tabix ${svcf}.nosr.snps.recode.vcf.gz
    #bcftools sort -T \$PWD ${svcf}.nosr.snps.recode.vcf.gz -Oz -o ${svcf}.nosr.snps.sorted.vcf.gz
    #bcftools index ${svcf}.nosr.snps.sorted.vcf.gz

    vcftools --vcf ${svcf}.notd.vcf \
      --remove-indels --recode --recode-INFO-all \
      --chr chr${chr} \
      --out ${svcf}.notd.snps
    sed -i "s/sample/${sample}_hap${hap}/" ${svcf}.notd.snps.recode.vcf
    bgzip ${svcf}.notd.snps.recode.vcf
    tabix ${svcf}.notd.snps.recode.vcf.gz
    bcftools sort -T \$PWD ${svcf}.notd.snps.recode.vcf.gz -Oz -o ${svcf}.notd.snps.sorted.vcf.gz
    bcftools index ${svcf}.notd.snps.sorted.vcf.gz

    """
}


process MERGE_VCF {
    fair true
    label 'merge_vcf'

    input:
    tuple val(chr), path(vcfs), path(vcf_index)

    output:
    tuple val(chr), path("all_sample_chr${chr}.syri.notd.snps.merged.vcf.txt")

    script:
    """
    #bcftools merge ${vcfs} -Oz -o all_sample_chr${chr}.syri.notd.snps.merged.vcf.gz
    # bcftools index all_sample_chr${chr}.syri.notd.snps.merged.vcf.gz
    # grep -v '^##' all_sample_chr${chr}.syri.notd.snps.merged.vcf.gz
    bcftools merge ${vcfs} \
    | grep -v '^##' \
    | cut --complement -f 3,8 \
    | grep -v ',' \
    > all_sample_chr${chr}.syri.notd.snps.merged.vcf.txt

    #awk 'NR>1 { print \$1"\t"\$2-1"\t"\$2}' all_sample_chr${chr}.syri.notd.snps.merged.vcf.txt \
    #> all_sample_chr${chr}.syri.notd.snps.merged.vcf.regions_for_tabix.txt
    """
}


process MERGE_VCF_TABIX {
    fair true
    label 'MERGE_VCF_TABIX'

    input:
    tuple val(chr), path(vcftxt)

    output:
    tuple val(chr), path("all_sample_chr${chr}.syri.notd.snps.merged.vcf.regions_for_tabix.txt")
//    path "all_sample_chr${chr}.syri.notd.snps.merged.vcf.gz.csi"

    script:
    """
    awk 'NR>1 { print \$1"\t"\$2-1"\t"\$2}' $vcftxt \
    > all_sample_chr${chr}.syri.notd.snps.merged.vcf.regions_for_tabix.txt
    """
}


process SYRI_ANNO {
    fair true
    label 'SYRI_ANNO'

    input:
    tuple val(chr), path(soutbed), path(soutbedtbi), path(mergedvcf)

    output:
    path "${soutbed.baseName}.notd.snp_anno.txt"

    script:
    """
    tabix ${soutbed} -R ${mergedvcf} \
    | sort -k1,1 -k2,2n \
    | uniq > ${soutbed.baseName}.notd.snp_anno.txt
    """
}

process DEL_MARKERS{
    fair true
    label 'SYRI_ANNO'

    input:
    tuple val(chr), path(soutbed), path(soutbedtbi), path(mergedvcf)

    output:
    path "${soutbed.baseName}.notd.snp_anno.txt"

    script:
    """

for c in {01..12}; do
for c in {02..12}; do
  cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr$c
  chrid=$c
  export chrid
  sbatch ../../tool/sbatch_files/SBATCH_get_deletion_markers.sh
done

cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/chr${chrid}
cd $cwd

srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=${SLURM_MEM_PER_CPU} \
python ../../tool/get_deletion_markers.py $cwd chr${chrid} ${SLURM_CPUS_PER_TASK}

from util import updatesnpgenotypes
import sys
"""
Required arguments:
pwd: Working directory
c: chromosome ID
nproc: Number of processes
"""
updatesnpgenotypes(*sys.argv[1:])


    echo 1

    """

}


workflow{
	a = CREATECHROMINDICES(references.concat(chromosomes))
	sams = MINIMAP2_CHR_TO_CHR(indexes)
	bams = SAM_TO_BAM(indexes, sams)
	(souts, svcfs) = RUN_SYRI(indexes, bams)
//	souts.view()
	sout_bed = SYRI_INDEX(souts)
//    sout_bed.view()
	svcf_filt = SYRI_FILTER(indexes, svcfs)
//	(merged_vcf, merged_vcf_index) = svcf_filt
	merged_vcf_txt = svcf_filt
			| groupTuple()
            | MERGE_VCF
    merged_vcf_tabix = MERGE_VCF_TABIX(merged_vcf_txt)

//    merged_vcf_txt.view()

    sout_bed_grouped = sout_bed
    | map { bed, tbi ->
      def key = bed.name.toString().tokenize('chr').get(1).tokenize('_').get(0)
      return tuple(key, bed, tbi)
    }

    sout_bed_grouped
    | combine(merged_vcf_tabix, by: 0)
    | set{sout_merged_vcf}

//    sout_merged_vcf.view()
    syri_anno =  SYRI_ANNO(sout_merged_vcf)
    syri_anno.view()









}



// USE working directory as $SCRATCH_DSS/pantohap
