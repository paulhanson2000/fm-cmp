library(data.table)

samp <- fread("1000GP_Phase3/1000GP_Phase3.sample")
thresh <- as.list(40/table(samp$GROUP)) # Calculate MAF threshold equivalent to MAC<40 for each ancestry group

leg_files <- paste0("1000GP_Phase3/1000GP_Phase3_chr",1:22,".legend.gz")
legs <- lapply(leg_files,fread) # Could add select=... but faster to read whole thing b/c then data.table uses mmap() I think
mapply(legs,1:22, FUN=function(leg,i) {
  inds <- leg[, which( (AFR > thresh$AFR & AFR < 1-thresh$AFR) |
                       (AMR > thresh$AMR & AMR < 1-thresh$AMR) |
                       (EAS > thresh$EAS & EAS < 1-thresh$EAS) |
                       (EUR > thresh$EUR & EUR < 1-thresh$EUR) |
                       (SAS > thresh$SAS & SAS < 1-thresh$SAS) ) ]
  writeLines(as.character(inds), paste0("data/mac_gt40_indices_chr",i,".txt"))
  # Make a file of the indices of the variants with MAC>40, will use to filter the .hap and .leg files.
})


dir.create("data")
# --- Filter .hap files ---
system(paste("parallel -j4",
  "'gzip -cd 1000GP_Phase3/1000GP_Phase3_chr{}.hap.gz |",
  " awk \"NR==FNR {{a[\\$1]; next}} FNR in a\" data/mac_gt40_indices_chr{}.txt - |",
  " gzip > data/chr{}_mac_gt40.hap.gz' ::: {1..22}"
))

# --- Filter .leg files ---
# Just print header
system(paste("parallel -j4",
  "'gzip -cd 1000GP_Phase3/1000GP_Phase3_chr{}.legend.gz |",
  " head -1 |",
  " gzip > data/chr{}_mac_gt40.leg.gz' ::: {1..22}"
))

# Print the rows whose index+1 are in the indices file. The +1 accounts for the header line.
system(paste("parallel -j4",
  "'gzip -cd 1000GP_Phase3/1000GP_Phase3_chr{}.legend.gz |",
  " awk \"NR==FNR {{a[\\$1]; next}} (FNR+1) in a\" data/mac_gt40_indices_chr{}.txt - |",
  " gzip >> data/chr{}_mac_gt40.leg.gz' ::: {1..22}"
))

# .hap files for the larger chromosomes don't fit in memory.
# I tried the vroom package, which can lazily load large files.
# Unfortunately, it was much slower than just using the ugly bash commands above. :(
# (Bash commands also work without overflowing memory because they stream files line-by-line.)
#library(vroom)
#cmd <- "gzip -cd ~/Downloads/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz | awk 'NR==FNR {a[$1]; next} FNR in a' indices.txt - > subset_chr22.hap.gz"
#hap <- vroom("~/Downloads/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz",col_types='c')
#hap <- hap[mac_gt40,]
