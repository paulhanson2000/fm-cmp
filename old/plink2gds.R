# See SNPRelate vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
# Specifically the "Format conversion from PLINK tex/binary files" section
# PLINK format 1000 Genomes data downloaded from: https://ctg.cncr.nl/software/magma

library(SNPRelate)
dir_base <- "../data/reference_panels/1000_genomes/"

eas_bed <- paste0(dir_base, "eas/g1000_eas.bed")
eas_fam <- paste0(dir_base, "eas/g1000_eas.fam")
eas_bim <- paste0(dir_base, "eas/g1000_eas.bim")

eur_bed <- paste0(dir_base, "eur/g1000_eur.bed")
eur_fam <- paste0(dir_base, "eur/g1000_eur.fam")
eur_bim <- paste0(dir_base, "eur/g1000_eur.bim")

sas_bed <- paste0(dir_base, "sas/g1000_sas.bed")
sas_fam <- paste0(dir_base, "sas/g1000_sas.fam")
sas_bim <- paste0(dir_base, "sas/g1000_sas.bim")

snpgdsBED2GDS(eas_bed, eas_fam, eas_bim, paste0(dir_base, "eas/g1000_eas.gds"))
snpgdsBED2GDS(eur_bed, eur_fam, eur_bim, paste0(dir_base, "eur/g1000_eur.gds"))
snpgdsBED2GDS(sas_bed, sas_fam, sas_bim, paste0(dir_base, "sas/g1000_sas.gds"))

#snpgdsSummary(paste0(dir_base, "eas/g1000_eas.gds"))
#snpgdsSummary(paste0(dir_base, "eur/g1000_eur.gds"))
#snpgdsSummary(paste0(dir_base, "sas/g1000_sas.gds"))

#(eas_genofile <- snpgdsOpen(paste0(dir_base,"eas/g1000_eas.gds")))
#(eur_genofile <- snpgdsOpen(paste0(dir_base,"eur/g1000_eur.gds")))
#(sas_genofile <- snpgdsOpen(paste0(dir_base,"sas/g1000_sas.gds")))

#tmp <- read.gdsn(index.gdsn(eas_genofile, path="sample.annot"))

