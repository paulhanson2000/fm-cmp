library(data.table)

frq <- fread("ld_output/EUR_LD_frq.frq")
bim <- fread("ld_output/EUR_LD_ref.bim")
ld_og <- readBin("ld_output/EUR_LD.ld.bin", what="numeric", n=nrow(frq)^2, size=4)
ld    <- matrix(ld_og, nrow=nrow(frq), ncol=nrow(frq))

isSymmetric(ld)
# Matrix is symmetric, good

ld2 <- as.data.table(ld)
all(sapply(seq(1,nrow(frq)), function(i) ld2[[i,i]]==1))
rm(ld2)
# All diagonal entries are 1, good

          # v equiv to ld_og
writeBin(as.vector(ld), "ld_rewritten.ld.bin", size=4)
ld_rewritten <- readBin("ld_rewritten.ld.bin", what="numeric", n=nrow(frq)^2, size=4)
identical(ld_og,ld_rewritten)
# TRUE, good
