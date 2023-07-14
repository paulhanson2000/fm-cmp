library(SeqArray)

seqBED2GDS("../data/ref/1kg/plink_format/eas/g1000_eas.bed",
           "../data/ref/1kg/plink_format/eas/g1000_eas.fam",
           "../data/ref/1kg/plink_format/eas/g1000_eas.bim",
           "magma_eas.gds")
meas <- seqOpen("magma_eas.gds", allow.duplicate=T)
seqSummary(meas)
