process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
    tuple val(name), val(sample), path("${name}.rmdup.sorted.bam"), path("${name}.rmdup.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads | samblaster |samtools view -Sb - | sambamba sort /dev/stdin -o ${name}.rmdup.sorted.bam
        samtools index ${name}.rmdup.sorted.bam ${name}.rmdup.sorted.bai
	"""
}
