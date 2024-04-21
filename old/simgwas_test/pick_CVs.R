library(data.table)

samp <- fread("data/1000GP_Phase3.sample")
leg <- fread("data/chr22_mac_gt40.leg.gz")
sapply(unique(samp$GROUP), function(g) writeLines(samp[GROUP==g,ID], paste0("sample_ids-",g,".txt")) )

# .sample file missing 2nd header line, and SEX should be numbers. See https://www.cog-genomics.org/plink/2.0/formats#sample
fwrite(rbind(list("0","D","D","D"),samp[,SEX:=ifelse(SEX=="male",1,2)]), "data/1000GP_Phase3.sample2", sep=' ')

# Calc 50kb window LD for each ancestry
sapply(unique(samp$GROUP), function(g) {
  samp_id_fnm <- paste0("sample_ids-",g,".txt")
  system(paste("./plink2",
    "--haps   data/chr22_mac_gt40.hap.gz",
    "--legend data/chr22_mac_gt40.leg.gz 22",
    "--sample data/1000GP_Phase3.sample2",
    "--keep", samp_id_fnm,
    "--mac 40 minor", # I know I've already filtered by the looser threshold of "at least one ancestry must have > 40", but for good comparisons it'd be nice if BOTH ancestries have MAC>40.
    "--r-phased",
    "--ld-window-kb 50",
    "--ld-window-r2 0",
    "--out", g
))})

# Variant LD averages in a 50kb window for each ancestry
ld_avgs <- sapply(unique(samp$GROUP), simplify=F, function(g) {
  ld <- fread(paste0(g,".vcor"))
  ld_avg <- ld[,mean(abs(PHASED_R)),by=ID_A]
})
ld_avgs <- Reduce(function(x,y) x[y,on="ID_A"], ld_avgs)
names(ld_avgs) <- c("id",unique(samp$GROUP))

# EUR/AFR variants with low,medium,high average LD
eur_l_ld_ids <- ld_avgs[EUR<0.1,          id]
eur_m_ld_ids <- ld_avgs[EUR>0.2 & EUR<0.4,id]
eur_h_ld_ids <- ld_avgs[EUR>0.5,          id]
afr_l_ld_ids <- ld_avgs[AFR<0.2,          id]
afr_m_ld_ids <- ld_avgs[AFR>0.2 & AFR<0.4,id]
afr_h_ld_ids <- ld_avgs[AFR>0.5,          id]

# I want to test three cases:
# 1. Variant w/ low LD in both ancestires (nearly independent)
# 2. Variant in high LD in both ancestries (hard to pick out, but the LD patterns across ancestries might still be different, even if they are both high)
# 3. Variant in higher LD in one ancestry than the o/ (ideal scenario for multi-ancestry fine-mapping)
# Variant should have roughly equal MAF in both ancestries s.t. it's really the LD that makes the difference.
v <- ld_avgs[id %in% intersect(eur_l_ld_ids, afr_l_ld_ids), .(id,AFR,EUR)] # rs117142056:39818793:A:G  LDavg: AFR:0.052 EUR:0.044;  MAF: AFR
v <- cbind(v, leg[id %in% v$id, .(AFR_maf=AFR,EUR_maf=EUR)])
# Options for 1.
# id,AFR,EUR,AFR_maf,EUR_maf,EUR-AFR,EUR_maf-AFR_maf
#* rs12169200:27638781:A:G,0.059,0.063,0.32,0.32,0.01,-0.00
# rs3076533:48883825:AGTGT:A,0.094,0.095,0.70,0.71,-0.02,-0.00

v <- ld_avgs[id %in% intersect(eur_h_ld_ids, afr_h_ld_ids), .(id,AFR,EUR)] # rs117142056:39818793:A:G  LDavg: AFR:0.052 EUR:0.044;  MAF: AFR
v <- cbind(v, leg[id %in% v$id, .(AFR_maf=AFR,EUR_maf=EUR)])
# Options for 2.
# id,AFR,EUR,AFR_maf,EUR_maf,EUR-AFR,EUR_maf-AFR_maf
#* rs5755304:35201115:T:C,0.51,0.53,0.66,0.61,0.03,-0.06
# rs5755301:35200936:T:G,0.52,0.56,0.66,0.61,0.04,-0.06

v <- ld_avgs[id %in% intersect(eur_h_ld_ids, afr_l_ld_ids), .(id,AFR,EUR)] # rs117142056:39818793:A:G  LDavg: AFR:0.052 EUR:0.044;  MAF: AFR
v <- cbind(v, leg[id %in% v$id, .(AFR_maf=AFR,EUR_maf=EUR)])
# Options for 3.
# id,AFR,EUR,AFR_maf,EUR_maf,AFR-EUR,AFR_maf-EUR_maf
# rs9608946:30892255:A:G,0.18,0.58,0.23,0.23,-0.40,-0.01
#* rs28625949:24655227:G:A,0.13,0.56,0.35,0.35,-0.44,-0.01



# --- OLD ---
#d <- fread("plink2.vcor")
#summary(abs(d$PHASED_R))
#avgs <- d[,mean(abs(PHASED_R)),by=ID_A]
#fwrite(avgs, "avg_lds-50kb_window-mac_gt40.csv.gz")
#png(); hist(avgs$V1); dev.off()

#sum(              avgs$V1>=0.5)/length(avgs$V1)
#sum(avgs$V1>0.1 & avgs$V1<0.5)/length(avgs$V1)
#sum(avgs$V1<=0.1            )/length(avgs$V1)
