Vcfs2Gds <- function(files, output_name, sample_ids=NULL, variant_ids=NULL, regions=NULL, bcftools_bin="bcftools") {
  #' Vcfs2Gds
  #'
  #' Convert a bunch of VCF/BCF files into one SeqArray GDS file.
  #'
  #' @param files VCF/BCF file path(s) or URL(s). Or, the path to a file listing one VCF/BCF filename or URL per line.
  #' @param output_name Filename of output GDS file. 
  #' @param sample_ids Character vector of sample IDs, to subset the VCFs by. Or, the path to a file listing one sample ID per line.
  #' @param variant_ids Character vector of variant IDs, to subset the VCFs by. Or, the path to a file listing one variant ID per line. TODO: only rs ids work currently.
  #' @param regions Character vector of regions e.g. "1:234-567". Or, the path to a UCSC BED format file. Or, a data.frame of bed-like format (1st column chrom, 2nd chromStart, 3rd chromEnd)
  #' @param bcftools_bin Path to the program "bcftools".
  #' @section A tip to increase speed: Even if you already provide variant_ids, also providing regions will speed things up a lot. It is easier for computers to compare chromosome and position numbers than ID strings.
  #' @section TODO: Support and chr:pos(:ref:alt) IDs? Support mixed types of IDs? Infer regions if all IDs are chr:pos(:ref:alt), and maybe mention in these docs that that optimization is made once implemented. Deal with "chr#" vs. "#" syntax for the user of not matching (I think bcftools norm or  bcftools annotate --rename-chrs can help). Allow GRanges regions input?

  if(system(bcftools_bin)==127) stop("Could not run ",bcftools_bin,". Please ensure bcftools is installed, or change the bcftools_bin argument for a typo.")
  if(!length(variant_ids)==1 || !file.exists(variant_ids)) { if(!all(grepl("^rs",variant_ids))) stop("<TODO> Only rs IDs are supported. In the meantime you could leave it blank (meaning all variants will be gotten), and then filter the resulting GDS file how you like.") }

  # If needed, convert the function inputs to files (makes the bcftools call simpler later)
  is_file <- function(x) length(x)==1 && file.exists(x)
  if(!is_file(      files) && !is.null(      files)) writeLines(      files, "/tmp/vcf_or_bcf_files.txt")
  if(!is_file( sample_ids) && !is.null( sample_ids)) writeLines( sample_ids, "/tmp/sample_ids.txt"      )
  if(!is_file(variant_ids) && !is.null(variant_ids)) writeLines(variant_ids, "/tmp/variant_ids.txt"     )
  if(!is_file(    regions) && !is.null(    regions)) {
    if(!is.data.frame(regions)) { # If user did not already input a bed-like data.frame
      a <- strsplit(regions, ":|-")
      regions <- data.frame(
        chrom      = lapply(a, '[', 1),
        chromStart = lapply(a, '[', 2),
        chromEnd   = lapply(a, '[', 3))
    }
    write.table(regions, "/tmp/regions.bed", row.names=F, col.names=F, quote=F)
  }

  # Creating command:
  # bcftools concat -Ou -f <files> -R <regions> | bcftools view -Ou -i <variant_ids> -S <sample_ids>
  bcftools_cmd <- paste(bcftools_bin, "concat -Ou -f /tmp/vcf_or_bcf_files.txt")
  if(!is.null(regions)) bcftools_cmd <- paste(bcftools_cmd, "-R /tmp/regions.bed")
  bcftools_cmd <- paste(bcftools_cmd,"|",bcftools_bin,"view -Ou")
  if(!is.null( variant_ids)) bcftools_cmd <- paste(bcftools_cmd, "-i ID=@/tmp/variant_ids.txt")
  if(!is.null(  sample_ids)) bcftools_cmd <- paste(bcftools_cmd, "-S /tmp/sample_ids.txt"     )
  if(!is.null(chr_name_map)) bcftools_cmd <- paste(bcftools_cmd, "|",bcftools_bin,"annotate --rename-chrs",chr_name_map)
  return(bcftools_cmd)
  
  vcf_connection <- pipe(bcftools_cmd)
  f <- seqVCF2GDS(vcf_connection, output_name)
  close(vcf_connection)
  unlink("/tmp/vcf_or_bcf_files.txt", "/tmp/sample_ids.txt", "/tmp/variant_ids.txt", "/tmp/regions.bed")
  f
}
