library(data.table)
library(SeqArray)
# See ?seqRecompress: VCF2GDS conversion takes tons of memory w/ default compression.
# Better to have lower compression during conversion, then use seqRecompress() after.
seqVCF2GDS("../data/ref/topmed/chr21.BRAVO_TOPMed_Freeze_8.vcf.gz", "../data/ref/topmed/topmed_chr21_freeze.gds", storage.option="ZIP_RA")
seqRecompress(                                                      "../data/ref/topmed/topmed_chr21_freeze.gds",                 "LZMA")
freeze21 <- seqOpen(                                                "../data/ref/topmed/topmed_chr21_freeze.gds")

# Oh, it can't be converted, it's not VCF I didn't notice :V
#seqVCF2GDS("../data/ref/topmed/chr21.BRAVO_TOPMed_coverage_hg38.txt.gz", "../data/ref/topmed/topmed_chr21_coverage.gds", storage.option="ZIP_RA")
#seqRecompress(                                                           "../data/ref/topmed/topmed_chr21_coverage.gds",                 "LZMA")
#freeze21 <- seqOpen(                                                     "../data/ref/topmed/topmed_chr21_coverage.gds")

tmp <- fread("../data/ref/topmed/chr21.BRAVO_TOPMed_coverage_hg38.txt", nrows=3)
