library(data.table)

frq <- fread("ld_output/EUR_LD_frq.frq")
bim <- fread("ld_output/EUR_LD_ref.bim")
ld <- readBin("ld_output/EUR_LD.ld.bin", what="numeric", n=nrow(frq), size=4)
