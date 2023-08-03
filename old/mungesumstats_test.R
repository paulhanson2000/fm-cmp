library(data.table)
library(MungeSumstats)

BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
reformatted <- MungeSumstats::format_sumstats(path=eduAttainOkbayPth, ref_genome="GRCh37")
example_sumstats <- fread(reformatted)

my_sumstats_path <- "fm-cmp/data/DIAMANTE2022/sumstat/DIAMANTE-EUR.sumstat.txt"
my_reformatted <- MungeSumstats::format_sumstats(path=my_sumstats_path, ref_genome="GRCh37")
# Started @ 12:54, failed 1:05

sumstat <- fread(my_sumstats_path)
setnames(sumstat, old=c("rsID", "chromosome(b37)", "position(b37)", "effect_allele", "other_allele", "effect_allele_frequency", "Fixed-effects_beta", "Fixed-effects_SE", "Fixed-effects_p-value"),
                  new=c("SNP",  "CHR",             "BP",            "A2",            "A1",           "FRQ",                     "BETA",               "SE",               "P"                    ))

my_reformatted <- MungeSumstats::format_sumstats(path=sumstat, ref_genome="GRCh37")
# started @ 1:18,
