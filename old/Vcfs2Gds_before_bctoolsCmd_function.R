Vcfs2Gds <- function(files, output_name, sample_ids=NULL, variant_ids=NULL, regions=NULL, chr_nm_map=NULL, info_fields=NULL, format_fields=NULL, bcftools_bin="bcftools") {
  #' Vcfs2Gds
  #'
  #' Convert a bunch of VCF files into one SeqArray GDS file.
  #'
  #' @param files VCF file path(s) or URL(s).
  #' @param output_name Filename of output GDS file. 
  #' @param sample_ids Character vector of sample IDs, to subset the VCFs by.
  #' @param variant_ids Character vector of variant IDs, to subset the VCFs by. TODO: only rs ids work currently.
  #' @param regions A data.frame of bed-like format (1st column chrom, 2nd chromStart, 3rd chromEnd, other cols unused)
  #' @param chr_nm_map A data.frame with two columns. Used to convert the VCF files' chromosome/contig names from the names in the first column to the second. E.g. to convert from "chr1" to "1".
  #' @param bcftools_bin Path to the program "bcftools".
  #' @section A tip to increase speed: Even if you already provide variant_ids, also providing regions will speed things up a lot. It is easier for computers to compare chromosome and position numbers than ID strings.
  #' @section TODO: Support and chr:pos(:ref:alt) IDs? Support mixed types of IDs? Infer regions if all IDs are chr:pos(:ref:alt), and maybe mention in these docs that that optimization is made once implemented. Deal with "chr#" vs. "#" syntax for the user of not matching (I think bcftools norm or  bcftools annotate --rename-chrs can help). Allow GRanges regions input? Support BCF files too?--Note would have to use pipe(..., "rb") maybe.

  if(system(bcftools_bin, ignore.stderr=T)==127) stop("Could not run ",bcftools_bin,". Please ensure bcftools is installed, or check the bcftools_bin argument for a typo.")
  if(!all(grepl("^rs",variant_ids))) stop("<TODO> Only rs IDs are supported. In the meantime you could leave it blank (meaning all variants will be gotten), and then filter the resulting GDS file how you like afterwards.")

                            writeLines(      files, "/tmp/vcf_files.txt"  )
  if(!is.null( sample_ids)) writeLines( sample_ids, "/tmp/sample_ids.txt" )
  if(!is.null(variant_ids)) writeLines(variant_ids, "/tmp/variant_ids.txt")
  if(!is.null( chr_nm_map)) write.table(chr_nm_map, "/tmp/chr_nm_map.txt", sep=' ', row.names=F, col.names=F, quote=F)
  if(!is.null(    regions)) {
    regions <- regions[order( unlist(regions[[1]]), unlist(regions[[2]]) )] # Sort by chrom,chromStart

    vcfs_contigs <- unique(c(sapply(files, function(f) {
      #       Gets   vvvvv   from VCF headers
      # ##contig=<ID=chr12,length=248956422,assembly=gnomAD_GRCh38>
      system(paste(bcftools_bin,"head",f,"| grep '##contig' | sed -e 's/.*ID=//' -e 's/,.*//'"), intern=T, ignore.stderr=T)
    })))

    if(!is.null(chr_nm_map)) {
      regions[[1]] <- sapply(regions[[1]], function(nm) {
        if(nm %in% vcfs_contigs)
          nm
        else if(nm %in% chr_nm_map[[2]])
          chr_nm_map[[1]][chr_nm_map[[2]]==nm] # Reverse-map regions' contig names to match the VCF
      })
    }

    if(is.null(chr_nm_map)) omitted_contigs <- vcfs_contigs[!(vcfs_contigs %in% c(regions[[1]]                    ))]
    else                    omitted_contigs <- vcfs_contigs[!(vcfs_contigs %in% c(regions[[1]], unlist(chr_nm_map)))]
    if(length(omitted_contigs) > 0)
      message("FYI, these contigs/chromosomes in the VCF files will be omitted because they are not mentioned in regions nor chr_nm_map:\n",
              paste(collapse='\n', omitted_contigs))

    write.table(regions, "/tmp/regions.bed", sep='\t', row.names=F, col.names=F, quote=F)
  }

  # Creating command:
  # bcftools concat -f <files> -R <regions>
  # -Ou | bcftools query -f TODO:
  # -Ou | bcftools view -i <variant_ids> -S <sample_ids>
  bcftools_cmd <- paste(bcftools_bin, "concat -f /tmp/vcf_files.txt")
  if(!is.null(    regions)) bcftools_cmd <- paste(bcftools_cmd, "-R /tmp/regions.bed -a")
  if(!is.null(
  bcftools_cmd <- paste(bcftools_cmd,"-Ou |",bcftools_bin,"view")
  if(!is.null(variant_ids)) bcftools_cmd <- paste(bcftools_cmd, "-i ID=@/tmp/variant_ids.txt")
  if(!is.null( sample_ids)) bcftools_cmd <- paste(bcftools_cmd, "-S /tmp/sample_ids.txt"     )
  if(!is.null( chr_nm_map)) bcftools_cmd <- paste(bcftools_cmd,"-Ou |",bcftools_bin,"annotate --rename-chrs /tmp/chr_nm_map.txt")
  bcftools_cmd <- paste(bcftools_cmd,"-o out.vcf")
  return(bcftools_cmd)
  system(bcftools_cmd)
  #bcftools_cmd <- paste(bcftools_cmd,"-Ov")
  
  vcf_connection <- pipe(bcftools_cmd, "rt")
  f <- seqVCF2GDS(vcf_connection, output_name)
  close(vcf_connection)
  unlink("/tmp/vcf_files.txt", "/tmp/sample_ids.txt", "/tmp/variant_ids.txt", "/tmp/regions.bed")
  f
}
