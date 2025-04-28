process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
      
        label "l_cpu"
        label "l_mem"
 
	
	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
    tuple val(name), val(sample), path("${name}.mdup.sorted.bam"), path("${name}.mdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex}.fa $fwd $rev | samblaster |samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.mdup.sorted.bam
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
        source activate gatk4610
        gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.ref}.fa -I $bam -L ${params.varbed}  --dont-use-soft-clipped-bases true -A StrandBiasBySample -minimum-mapping-quality 0 --mapping-quality-threshold-for-genotyping 0 --enable-dynamic-read-disqualification-for-genotyping true --flow-filter-alleles-qual-threshold 0 -O ${name}.vcf
        """
}

process NORMALIZACE {
        tag "NORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(gatk)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo NORMALIZACE $name
        bcftools norm -m-both -f ${params.ref}.fa -o ${name}.norm.vcf $gatk
        bgzip ${name}.norm.vcf
        tabix ${name}.norm.vcf.gz
        """
}

process ANOTACE_ACGT {
        tag "ANOTACEACGT on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.acgt.vcf")

        script:
        """
        source activate gatk4610
        echo ANOTACEACGT $name
        gatk --java-options "-Xmx4g"  VariantAnnotator   -V ${name}.norm.vcf.gz -O ${name}.norm.acgt.vcf.gz --resource:ACGT ${params.ACGT} --expression ACGT.AF   --expression ACGT.AC   --expression ACGT.AC_Hom   --expression ACGT.AC_Het   --expression ACGT.AC_Hemi

        gunzip ${name}.norm.acgt.vcf.gz
        """
}

process ANOTACE_OMIM {
        tag "ANOTACEOMIM on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(acgt)

        output:
        tuple val(name), val(sample), path("${name}.norm.vaf.omim.vcf")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACEOMIM $name
        bcftools norm -m-both -f ${params.ref}.fa -o ${name}.norm.vcf $acgt

        bcftools view ${name}.norm.vcf  -o ${name}.pom.bcf
        bcftools +fill-tags ${name}.pom.bcf -Ob -o ${name}.pom2.bcf -- -t FORMAT/VAF
        bcftools convert -O v -o ${name}.norm.vaf.vcf ${name}.pom2.bcf

bcftools annotate -a ${params.AR} -h ${params.ARheader} -c CHROM,FROM,TO,dedicnostAR -l dedicnostAR:append -m -xx ${name}.norm.vaf.vcf > ${name}.pom3
bcftools annotate -a ${params.AD} -h ${params.ADheader} -c CHROM,FROM,TO,dedicnostAD -l dedicnostAD:append -m -aa ${name}.pom3 > ${name}.pom4
bcftools annotate -a ${params.Xlinked}  -h ${params.Xlinkedheader} -c CHROM,FROM,TO,dedicnostXlinked -l dedicnostXlinked:append -m -bb ${name}.pom4 > ${name}.pom5
bcftools annotate -a ${params.Ylinked}  -h ${params.Ylinkedheader} -c CHROM,FROM,TO,dedicnostYlinked -l dedicnostYlinked:append -m -cc ${name}.pom5 > ${name}.pom6

bcftools annotate -a ${params.omimphenotyp} -h ${params.omimphenotypheader} -c CHROM,FROM,TO,fenotyp -l fenotyp:append -m -yy ${name}.pom6 | sed  's/xx/dedicnostAR=0/' | grep -v 'INFO=<ID=dedicnostAR=0' | sed  's/yy/fenotyp=0/' | grep -v 'INFO=<ID=fenotyp=0' | sed  's/aa/dedicnostAD=0/' | grep -v 'INFO=<ID=dedicnostAD=0' | sed  's/bb/dedicnostXlinked=0/' | grep -v 'INFO=<ID=dedicnostXlinked=0' | sed  's/cc/dedicnostYlinked=0/' | grep -v 'INFO=<ID=dedicnostYlinked=0' > ${name}.norm.vaf.omim.vcf

        """
}

process ANOTACE_annovar {
       tag "ANOTACE on $name"
       //publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(normalizovany)

        output:
        tuple val(name), val(sample), path("${name}.norm.vaf.omim.vcf.hg38_multianno.vcf.gz"), path("${name}.norm.vaf.omim.vcf.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput $normalizovany ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20240917,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 20 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.vaf.omim.vcf.hg38_multianno.vcf
        tabix ${name}.norm.vaf.omim.vcf.hg38_multianno.vcf.gz
        """
}


process VCF2TXT {
       tag "VCF2TXT on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
       
        input:
        tuple val(name), val(sample), path("${name}.norm.vaf.omim.vcf.hg38_multianno.vcf.gz"), path("${name}.norm.vaf.omim.vcf.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.txt")

        script:
        """
        echo VCF2TXT $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}.fa  --show-filtered  -V ${name}.norm.vaf.omim.vcf.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -GF VAF -F dedicnostAR -F dedicnostAD -F dedicnostXlinked -F dedicnostYlinked -F fenotyp  -F ACGT.AF -F ACGT.AC -F ACGT.AC_Hom -F ACGT.AC_Het -F ACGT.AC_Hemi -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.txt
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

