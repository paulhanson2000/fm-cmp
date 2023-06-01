library(data.table)
library(Rsamtools)
library(bedr)

tbx <- TabixFile("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
headerTabix(tbx)
countTabix(tbx, param=GRanges(c("12"), IRanges(start=c(27630940), end=c(27830940))) )
reg <- scanTabix(tbx, param=GRanges(c("12"), IRanges(start=c(27630940), end=c(27830940))) )
reg <- fread(tteaxtConnection(reg), sep='\t', header=F)
