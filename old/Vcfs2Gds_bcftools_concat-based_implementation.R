Vcfs2Gds <- function(files, sample_ids, variant_ids, regions, output_name, n_thread=parallel::detectCores(), scratch_dir="/tmp/") {
  if(!length(variant_ids)==1 || !file.exists(variant_ids)) { if(!all(grepl("^rs",variant_ids))) stop("Only rs IDs are supported. TODO") }

  # If needed, convert the function inputs to files (makes the bcftools call simpler later)
  if(!length(      files)==1 || !file.exists(      files)) writeLiens(      files,   "/tmp/vcf_or_bcf_files.txt")
  if(!length( sample_ids)==1 || !file.exists( sample_ids)) writeLines( sample_ids,  "/tmp/sample_ids.txt")
  if(!length(variant_ids)==1 || !file.exists(variant_ids)) writeLines(variant_ids, "/tmp/variant_ids.txt")
  if(!length(    regions)==1 || !file.exists(    regions)) {
    if(!is.data.frame(regions)) {
      a <- strsplit(regions, ":|-")
      bed <- data.frame(
        chrom      = lapply(a, '[', 1),
        chromStart = lapply(a, '[', 2),
        chromEnd   = lapply(a, '[', 3))
      write.table(bed, "/tmp/regions.bed")
  }


  # Faster / more readable to do filtering by variant id after converted to seqArray?
  bcftools_cmd <- paste(
    "bcftools concat",
    "-f /tmp/vcf_or_bcf_files.txt",
    "-R /tmp/regions.bed",
    "-Ou", # (Optimization) Output uncompressed here since piping to another bcftools call
    "--naive", # TODO: only add this conditionally
    "|",
    "bcftools view",
    "-o",output_name,
    "-S /tmp/sample_ids.txt",
    "-i ID=@/tmp/variant_ids.txt" # TODO: only add this conditionally
    "--threads ", n_thread, # Don't need threads in the concat part b/c bcftools only multithreads file compression.
  )

  # It would be more memory-efficient to call bcftools once per file and then seqVCF2GDS it immediately before going to the next.
  # Then at the end, seqMerge all of them.
    # I was afraid it would be less readable, but actually if it would remove the piping nonsense then it'd be fine.
    # Although it would require keeping track of a bunch of indiv gds files to merge at the end, requiring the additiong of a new scratch_dir param...

  # TODO: oh right ano/ thing that was bothering me was that this function should also accept GDS file right, actually, well, should it? Maybe not.
    # Maybe this function's responsibility is not to "subset user's input no matter what form it's in and output a GDS", but instead more simply:
      # "Take VCF/BCFs and concat and output a GDS w/ optional filtering (which is recommemnded to do now b/c bcftools streams)"
    

  system(bcftools_cmd)
  # TODO: if fail, try again removing --naive


    
}
