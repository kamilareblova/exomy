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

23.4. 2026 

pridany sloupec z vep anotace VEP_consequence

-----------------------------------------------------------------------------

29.4. - update omim databaze - nove stazeno + pridani sloupce Mane select z Vep anotace

--------------------------------------------------------------

12.5. - pridana analyza LOH a GNUPLOT

----------------

15.5. - opraveni grep za grep -w u lohu

--------------------------------------------

19.6. - pridani qc kontroly pro mapping pomoci picardu.
vysledky z qualimap do samostatneho adresare qualimap
a z obou adresaru 2 multiqc reporty

------------------------------------------------

22.6. - update databaze v splice ai na spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz
zmena v nextflow.config
-------------------------------------------------
