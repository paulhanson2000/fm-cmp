library(SeqArray)
gnomad_chr12_url <- "gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr12.vcf.bgz"
bcftools_cmd <- paste("bcftools view -r chr12:30764085-30774085 -o my.vcf", gnomad_chr12_url)
system(bcftools_cmd)

seqVCF2GDS("my.vcf", "out.gds")

bcftools_cmd <- paste("bcftools view -r chr12:30764085-30774085", gnomad_chr12_url)
vcf_con <- pipe(bcftools_cmd, "rt")
seqVCF2GDS(vcf_con, "out.gds")
