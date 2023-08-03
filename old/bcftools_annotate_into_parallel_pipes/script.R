{
files <- readLines("/tmp/vcf_files.txt")
pipes <- paste0("anno",1:length(files),".fifo")
writeLines(pipes,"my_pipes.txt")

mapply(files, pipes, 1:length(files), FUN=function(fnm,pnm,i) {
  system(paste("mkfifo",pnm))
  p <- parallel:::mcfork()
  if(inherits(p,"masterProcess")) {
    system(paste("bcftools annotate -a https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz -c ID",fnm,">",pnm))
  }
})

system("bcftools concat -f my_pipes.txt | bcftools head -n 10 > annotated.vcf")
}

# Doesn't work, bcftools annotate -a <dbSNP> -c ID causes file to need to be re-indexed I guess.
# I think it's for the best I can't get this to work lol.
