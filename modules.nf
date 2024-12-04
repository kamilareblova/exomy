process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
      
        label "l_cpu"
        label "l_mem"
 
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
    tuple val(name), val(sample), path("${name}.mdup.sorted.bam"), path("${name}.mdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex}.fa $reads | samblaster |samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.mdup.sorted.bam
        samtools index ${name}.mdup.sorted.bam ${name}.mdup.sorted.bai
	"""
}

process GATK {
       tag "GATK on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.vcf")

        script:
        """
        echo GATK $name
        source activate gatk4
        gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.ref}.fa -I $bam -L ${params.varbed}  --dont-use-soft-clipped-bases true -A StrandBiasBySample -minimum-mapping-quality 0 --mapping-quality-threshold-for-genotyping 0 --enable-dynamic-read-disqualification-for-genotyping true --flow-filter-alleles-qual-threshold 0 -O ${name}.vcf
        """
}

process NORMALIZACE {
        tag "NORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(gatk)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf")

        script:
        """
        source activate bcftoolsbgziptabix
        echo NORMALIZACE $name
        bcftools norm -m-both -f ${params.ref}.fa -o ${name}.norm.vcf ${gatk}.vcf
        """
}
process ANOTACE {
       tag "ANOTACE on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(normalizovany)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.hg38_multianno.vcf.gz"), path("${name}.norm.vcf.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput $normalizovany ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20221231,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 5 -exonicsplicing',,,,,,,,,,, --remove
        bgzip ${name}.norm.vcf.hg38_multianno.vcf
        tabix ${name}.norm.vcf.hg38_multianno.vcf.gz
        """
}

process VCF2TXT {
       tag "VCF2TXT on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.hg38_multianno.vcf.gz"), path("${name}.norm.vcf.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.txt")

        script:
        """
        echo VCF2TXT $name
        source activate java
        ${params.gatk} --java-options "-Xmx4g" VariantsToTable -R ${params.ref}.fa  --show-filtered  -V ${name}.norm.vcf.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.SB.txt
        """
}

process COVERAGE {
          tag "COVERAGE on $name"
       publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
         container "staphb/samtools:1.20"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.coveragefin.txt")

        script:
        """
        echo ANOTACE $name
        samtools bedcov ${params.varbed} $bam -d 20 > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print (\$6/(\$3-\$2))*100"%"}' ${name}.COV > ${name}-procento-nad-20
        paste ${name}.COV-mean ${name}-procento-nad-20 > vysledek
        echo "chr" "start" "stop" "name" ${name}.COV-mean ${name}-procento-nad-20 > hlavicka
        sed -i 's/ /\t/'g hlavicka
        paste ${params.varbed} vysledek > coverage
        cat hlavicka coverage > ${name}.coveragefin.txt
        sed -i -e "s/\r//g" ${name}.coveragefin.txt
        """
}

