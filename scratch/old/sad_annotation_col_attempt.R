# TODO: how to handle stuff like TFBSs where each singular row of the bed file oughta be its own annot and have its name in the bed file considered as well?
# TODO: bad code, loads all beds at once into memory, and they could be big and there could be a lot. Do one at a time.
appendAnnotCol <- function(sumstat, bed_file, annot_name) {
  bed <- fread(bed_file, select=c(1:3), col.names=c("chrom","chromStart","chromEnd"), sep=' ')
  bed <- bed[chrom %in% paste0("chr",1:22)] # TODO: temporary, only doing autosomes for now; and related: TODO: use factors for chr
  bed[, chrom := as.integer(sub("chr","",chrom))] # TODO: use factor for chr

  # Optimization (?) TODO: test if actually helps
  bed <- bed[ chrom %in% unique(sumstat$chr) &
              chromStart >= min(sumstat$pos) &
              chromEnd   <= max(sumstat$pos) ]

  #apply(bed, 1, function(r) {
  #  sumstat[ chr == r["chrom"]      &
  #           pos >= r["chromStart"] &
  #           pos <= r["chromEnd"]   ,
  #           (annot_name) := 1      ]
  #})

  annot_col <- c()
  for(v in 1:nrow(sumstat)) {
    vchr <- sumstat[v,chr]
    vpos <- sumstat[v,pos]
    chr_yes <- bed$chrom
    lpos_yes <- bed$chromStart <= vpos
    rpos_yes <- bed$chromEnd >= vpos
    yes <- any(mapply(chr_yes, lpos_yes, rpos_yes, FUN=all))
    annot_col <- c(annot_col,yes)
    #annot_col <- c(annot_col, c(bed$chrom      == vchr,
    #                                    bed$chromStart <= vpos,
    #                                    bed$chromEnd   >= vpos))
  }

  #mapply(sumstat$chr, sumstat$pos, FUN = ...
  
  annot_col
}
#appendAnnotCols <- Vectorize(

appendAnnotCols <- function(sumstat, bed_files) {
  beds <- lapply(bed_files, function(f) {
    bed <- fread(f, select=c(1:3), col.names=c("chrom","chromStart","chromEnd"), sep=' ') # bed files are standard, first three cols will always be this.
    bed <- bed[chrom %ni% c("chrX","chrY","chrM")] # TODO: temporary, only doing autosomes for now; and related: TODO: use factors for chr
    bed[, chrom := as.integer(sub("chr","",chrom))] # TODO: use factor for chr
  })
  sumstat <- cbind(sumstat, setnames(data.table(t(rep(0,length(beds)))),names(bed_files))) # Just adds 0-filled cols, one per bed file.
  
  mapply(beds, 1:length(beds), FUN = function(bed, bed_i) { 
    apply(bed, 1, function(r) {
      sumstat[ chr == r["chrom"]      & 
               pos >= r["chromStart"] &
               pos <= r["chromEnd"]   ,
               (bed_i) := 1           ]
  })})
  sumstat
}

# TODO: warning, some bed files may be ' ' OR '\t' delim'd, need s/t to guess which is the correct.
tmp_beds <- c("../data/annot/paintor_annot_lib/FANTOM5/cardiac_fibroblast_differentially_expressed_enhancers.bed",
              "../data/annot/paintor_annot_lib/Roadmap_ChromeHMM_15state/E001_15_coreMarks_mnemonics.bed.10_TssBiv.ES-I3_Cell_Line",
              "../data/annot/paintor_annot_lib/GeneElements_Gencode/GenCode.CDS.hg19")
names(tmp_beds) <- c("annotA","annotB","annotC")
#tmp <- appendAnnotCols(sumstats[[1]], tmp_beds)
tmp <- appendAnnotCol(sumstats[[1]], tmp_beds[3], "annotC")


if(!dir.exists("in/fgwas")) dir.create("in/fgwas")
if(!dir.exists("out/fgwas")) dir.create("out/fgwas")

mapply(sumstats, names(sumstats), FUN=function(s, nm) {
  s <- s[,    .(  rsid,  chr,  pos,  z,  se)]
  setnames(s, c("SNPID","CHR","POS","Z","SE")) # fgwas requires specific header names.
  s <- s[order(POS)] # fgwas requires ordering by pos
  s[, chr := paste0("chr",CHR)] # fgwas uses hg19 "chr#" syntax, not just the number.
  s[, `:=`(F=0, N=0)] # fgwas errors run if missing F & N cols, even though they're overridden by SE.

  filename <- paste0("in/fgwas/",nm,"-fgwas.gz")
  fwrite(s, filename, sep=' ', compress="gzip")
  system(paste("../third_party/fgwas/src/fgwas -i", filename,
               "-bed ../data/annot/paintor_annot_lib/FANTOM5/acinar_cell_differentially_expressed_enhancers.bed",
               "-o out/fgwas/tmp_test")) # TODO: TMP hardcoded test
  # TODO: problem w/ bed file seems to be the restriction that (Section 4.4): "all SNPs must fall within a segment defined in the file"
    # So it's prolly better just to add the annotation cols myself....
  # While you code this, think about the best interface for ano/ user adding annots would be. I guess kinda as it already is now; ask the user to supply bed files and the annots will be automatically added as cols. And maybe ALSO ask them if they already have annot cols in their substats? IDK, maybe not that. But if they DO have it already in there sumstats it would save uneccessary computation of adding them... idk.
})
