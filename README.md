# exomy
14.4. 2026 jsem pridala fastqc a multiqc report k workflow,
zaroven, se bude ukladat jen vcf, normalizvane vcf+index a merged file.
Flagstat jsem zakomentovala, pouziva se jen qualimap.

------------------------------------------------------------------------------
16.4. 2026

menim:
bcftools norm -m-both -f ${params.ref}.fa -o ${name}.norm.vcf ${name}.vaf.vcf

na

bcftools norm -f ${params.ref}.fa -m -both ${name}.vaf.vcf -o ${name}.norm.vcf

asi to nema vliv, porovnavala jsem rozdily, ale konvence je druhy zapis

------------------------------------------------------------------------------


22.4. 2026

pridane virtualni panely:
BRONCO  genoderm  HBOC  hemato  lipanek  LYNCH  pojiva  prekonc  skelet
------------------------------------------------------------------------------

