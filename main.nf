include { ALIGN; GATK; NORMALIZACE; ANOTACE_ACGT; ANOTACE_OMIM; ANOTACE_annovar; VCF2TXT; COVERAGE } from "${params.projectDirectory}/modules"

workflow {
        rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    .splitCsv(header: true)
    .map { row ->
        def baseDir = new File("${params.baseDir}")
        def runDir = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(row.run)
            }
        })[0] //get the real folderName that has prepended date

        def fileR1 = file("${runDir}/processed_fastq/${row.name}_R1.fastq.gz", checkIfExists: true)
        def fileR2 = file("${runDir}/processed_fastq/${row.name}_R2.fastq.gz", checkIfExists: true)

                def meta = [name: row.name, run: row.run]
        [
            meta.name,
            meta,
            fileR1,
            fileR2,
                ]
    }
     . view()

aligned = ALIGN(rawfastq)
varcalling = GATK(aligned)
normalizovany = NORMALIZACE(varcalling)
anotovanyacgt = ANOTACE_ACGT(normalizovany)
anotovanyomim = ANOTACE_OMIM(anotovanyacgt)
anotovany = ANOTACE_annovar(anotovanyomim)
anotovanyfin = VCF2TXT(anotovany)
coverage = COVERAGE(aligned)
}
