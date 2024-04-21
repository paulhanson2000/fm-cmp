library(qqman)
library(data.table)
setDTthreads(1)

gwas <- fread("plink.assoc")
gwas <- gwas[!is.na(P)]
highlight_rsid <- gwas[BP==2190692,SNP]

png("manhattan.png", width=800,height=400)
manhattan(gwas, highlight=highlight_rsid)
dev.off()

png("qq.png", width=400,height=400)
qq(gwas$P)
dev.off()
