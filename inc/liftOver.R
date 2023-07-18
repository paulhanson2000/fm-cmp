# This function can accept a file as input, and can be run standalone on the command line:
  # R -q -e "source(liftOver.R); liftOver("path/to/liftOver", "my_file.bed", "hg19", "hg38")
liftOver <- function(liftover_bin, to_lift, from_build, to_build, out_file=NULL) {
  # Get chain file
  chain_file_path <- dirname(liftover_bin)
  chain_file_name <- paste0(from_build, "To", tools::toTitleCase(as.character(to_build)), ".over.chain.gz")
  chain_file <- paste0(chain_file_path,"/",chain_file_name)

  if(from_build == to_build) stop("Using liftOver to convert from one reference genome to the same reference genome!")
  if(!file.exists(chain_file)) {
    err_code <- download.file(paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",from_build,"/liftOver/",chain_file_name), destfile = chain_file)
    if(err_code > 0) stop("Please download the chain file https://hgdownload.soe.ucsc.edu/goldenPath/",from_build,"/liftOver/",chain_file_name," and add it to ",chain_file_path,"/.")
  }

  wrote_tmp_file <- F
  if(!is.character(to_lift)) { # User passed in a data.frame
    wrote_tmp_file <- T
    write.table(to_lift, "/tmp/to_lift.bed", sep=' ', row.names=F, col.names=F, quote=F)
    to_lift <- "/tmp/to_lift.bed"
  }
  if(!file.exists(to_lift)) stop("File ",to_lift," does not exist!")

  if(is.null(out_file)) out_file <- "/tmp/lifted.bed"

  system(paste(liftover_bin, "-multiple -noSerial", to_lift, chain_file, out_file, "/tmp/unmapped.tmp"))

  n_bed_rows <- as.integer(system(paste("wc -l", to_lift, "| awk '{print $1}'"), intern=T))
  n_unmapped <- tryCatch(nrow(read.table("/tmp/unmapped.tmp")), error=function(e) 0)
  if(n_unmapped > 0) message(n_unmapped,"/",n_bed_rows," regions in ",to_lift," had parts which could not be liftOver'd:")
  system("cat /tmp/unmapped.tmp")

  out <- read.table(out_file)
  if(wrote_tmp_file) unlink(to_lift)
  unlink(c("/tmp/lifted.bed", "/tmp/unmapped.tmp"))
  out
}
