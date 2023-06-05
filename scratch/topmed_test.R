library(SeqArray)
seqVCF2GDS("../data/ref/topmed/chr21.BRAVO_TOPMed_Freeze_8.vcf.gz", "topmed_chr21_test.gds")
gds <- seqOpen("topmed_chr21_test.gds", allow.duplicate=T)
