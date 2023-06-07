# IDK kinda useless but maybe good reference
# And frankly it was a lot of typing so I don't want to throw it away if there's any chance I'll regret it

subset_sumstats <- list(list())
for(l  in 1:nrow(loci_configs)) {

  lchr     <- loci_configs[l, chr    ]
  pos_min <- loci_configs[l, pos_min]
  pos_max <- loci_configs[l, pos_max]

  for(ss in 1:nrow(data_configs)) {

    subset_sumstats[[l]][[ss]] <- sumstats[[ss]][chr %in% lchr &
                                                 pos > pos_min &
                                                 pos < pos_max ]
}} rm(l,ss)

# Or

subsetSumstats <- function(sumstats, chrs, pos_min, pos_max) {
  sumstats <- lapply(sumstats, '[', chr %in% chrs &
                                    pos > pos_min &
                                    pos < pos_max )
}
subsetSumstatsV <- Vectorize(subsetSumstats, vectorize.args=c("chrs","pos_min","pos_max"))
subset_sumstats <- subsetSumstatsV(sumstats, loci_configs$chr, loci_configs$pos_min, loci_configs$pos_max)

# Or

subset_sumstats <- lapply(1:nrow(locus_configs), function(r) {
                   lapply(         sumstats,     function(ss) {

  ss[chr == locus_configs[r,chr]     &
     pos >  locus_configs[r,pos_min] &
     pos <  locus_configs[r,pos_max] ]
})})
