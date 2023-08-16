bcftoolsPipableCmd <- function(files,
                               output_name=NULL, output_type=NULL,
                               regions=NULL, variant_ids=NULL, sample_ids=NULL,
                               query=NULL,
                               chr_nm_map=NULL, exclude_annos=NULL, extra_cmds=NULL,
                               scratch_dir="/tmp", bcftools_bin="bcftools") {
  #' bcftoolsPipableCmd
  #'
  #' A wrapper for a couple parts of bcftools. Combine and subset VCF/BCF files.
  #'
  #' @param files Character vector of VCF/BCF file paths and/or URLs. If multiple files, they will be concatenated, but must have the same sample columns. Combination will be much faster if no regions are given to subset by, the files are either all VCF or all BCF, and all their headers are compatible (useful to combine VCFs split by chromosome).
  #' @param output_name Filename of output file. Filename suffix determines the type (.vcf, .vcf.bgz, .bcf, .bcf.bgz).
  #' @param output_type "z"/"v" = compressed/uncompressed VCF, "b"/"u" = compressed/uncompressed BCF. Note that the file extension of output_name will take precedence.
  #' @param regions data.frame of bed-like format (1st column chromosome, 2nd start position (0-indexed), 3rd end, other cols unused)
  #' @param variant_ids Character vector of variant IDs, to subset the VCF/BCF(s) by.
  #' @param sample_ids Character vector of sample IDs, to subset the VCF/BCF(s) by.
  #' @param query <TODO document me>
  #' @param chr_nm_map data.frame with two columns. Used to convert the VCF files' chromosome/contig names from the names in the first column to the second. E.g. to convert from "chr1" to "1". You can assume this name mapping is applied before the VCF is subset by regions and/or variant_ids.
  #' @param info_fields <TODO> Character vector of INFO fields to keep, e.g. c("AC","AF","AN"). If left NULL, all fields are gotten.
  #' @param format_fields <TODO> Character vector of FORMAT fields to keep, e.g. c("GT","DP"). If left NULL, all fields are gotten.
  #' @param exclude_annos <Document TODO, maybe kinda redundant w/ info_fields and format_fields options? see docs of bcftools anno -x>
  #' @param extra_cmds <Document TODO, I know you've been thinking about rm'ing this but I think it's useful, since o/w would have to write to an intermediate file; the advantage of tacking on cmds like this is you can create one big single ~pipable~ cmd>
  #' @param scratch_dir Directory where temporary files should be written.
  #' @param bcftools_bin Path to the program "bcftools".
  #' @section A tip to increase speed: Even if you already provide variant_ids, also providing regions will speed things up a lot. It is easier for computers to compare chromosome and position numbers than ID strings.
  #' @section TODO: Support and chr:pos(:ref:alt) IDs? Support mixed types of IDs? Infer regions if all IDs are chr:pos(:ref:alt), and maybe mention in these docs that that optimization is made once implemented. Allow GRanges regions input?

  if(127 == system(bcftools_bin, ignore.stderr=T)) stop("Could not run ",bcftools_bin,". Please ensure bcftools is installed, or check the bcftools_bin argument for a typo.")
  if(!all(grepl("^rs",variant_ids))) stop("<TODO> Only rs IDs are supported. In the meantime you could leave it blank (meaning all variants will be gotten), and then filter the resulting GDS file how you like afterwards.")

  files_fnm       <- paste0(scratch_dir,"/files_list.txt")
  regions_fnm     <- paste0(scratch_dir,"/regions.bed")
  variant_ids_fnm <- paste0(scratch_dir,"/variant_ids.txt")
  sample_ids_fnm  <- paste0(scratch_dir,"/sample_ids.txt")
  chr_nm_map_fnm  <- paste0(scratch_dir,"/chr_nm_map.txt")

                            writeLines(      files,       files_fnm)
  if(!is.null(variant_ids)) writeLines(variant_ids, variant_ids_fnm)
  if(!is.null( sample_ids)) writeLines( sample_ids,  sample_ids_fnm)
  if(!is.null( chr_nm_map)) write.table(chr_nm_map,  chr_nm_map_fnm, sep=' ', row.names=F, col.names=F, quote=F)
  if(!is.null(    regions)) {
    #regions <- regions[order( unlist(regions[[1]]), unlist(regions[[2]]) )] # Sort by chrom,chromStart TODO: not necessary I don't think

    #       Gets   vvvvv   from VCF headers
    # ##contig=<ID=chr12,length=248956422,assembly=gnomAD_GRCh38>
    vcfs_contigs <- unique(unlist(sapply(files, function(f) {
      system(paste(bcftools_bin,"head",f,"| grep '##contig' | sed -e 's/.*ID=//' -e 's/>.*//' -e 's/,.*//'"), intern=T, ignore.stderr=T)
    })))

    # Use chr_nm_map to change regions contig names to match VCF's contig naming. Why not change the VCF's names first instead? See notes to self at the bottom.
    if(!is.null(chr_nm_map)) {
      regions[[1]] <- sapply(regions[[1]], function(nm) {
        # Map regions to preferred naming style. Then reverse-map to whatever is in VCF file.
        if(nm %in% unique(chr_nm_map[[1]]))
          nm <- chr_nm_map[[2]][chr_nm_map[[1]] == nm]

        if(!(nm %in% vcfs_contigs))  {# If not already matching VCF, then try reverse-mapping 
          nm2 <- chr_nm_map[[1]][chr_nm_map[[2]] %in% nm &
                                 chr_nm_map[[1]] %in% vcfs_contigs]
          if(length(nm2)==1) return(nm2) # i.e. length not 0
        }
        return(nm) # <TODO: better wording> Either nm is matching a contig name in VCF, or if not just return nm. If it doesn't match anything the region will just get dropped.
      })
    }

    if(is.null(chr_nm_map)) omitted_contigs <- vcfs_contigs[!(vcfs_contigs %in% c(regions[[1]]                    ))]
    else                    omitted_contigs <- vcfs_contigs[!(vcfs_contigs %in% c(regions[[1]], unlist(chr_nm_map)))]
    if(length(omitted_contigs) > 0)
      message("FYI, these contigs/chromosomes in the VCF files will be omitted because they are not mentioned in regions nor chr_nm_map:\n",
              paste(collapse='\n', omitted_contigs))

    write.table(regions, regions_fnm, sep='\t', row.names=F, col.names=F, quote=F)
  }


  # Now to paste() together a command. Something like:
  # bcftools concat -f <files> -R <regions>
  # -Ou | bcftools view -i <variant_ids> -S <sample_ids>
  # -Ou | bcftools annotate --rename-chrs <chr_nm_map>

  # Shorthand to avoid typing "bcftools_cmd <- paste(bcftools_cmd, thing_to_append)"
  bcftools_cmd <- ""
  a <- function(...) bcftools_cmd <<- paste0(bcftools_cmd, ...)

  # bcftools concat
  a(bcftools_bin, " concat -f ",files_fnm)
  if(!is.null(regions))
    a(" -a -R ",regions_fnm)

  # TODO: possible performance improvement using "--naive" concat option, but causing more issues than it's worth atm
  #headers_compatible <- 0 == system(paste(bcftools_bin,"concat --naive -f",files_fnm,"| bcftools head"), ignore.stderr=T, ignore.stdout=T)
  #if(is.null(regions) && headers_compatible)
  #  a(" --naive ")
  
  # bcftools view
  a(" -Ou | ",bcftools_bin," view")
  if(!is.null(variant_ids))
    a(" -i ID=",variant_ids_fnm)
  if(!is.null(sample_ids))
    a(" -S ",sample_ids_fnm)

  # bcftools annotate
  a(" -Ou | ",bcftools_bin," annotate")
  if(!is.null(chr_nm_map))
    a(" --rename-chrs ",chr_nm_map_fnm)
  if(!is.null(exclude_annos))
    a(" -x '", paste(collapse=',',exclude_annos), "'")

  #bcftools query
  if(!is.null(query))
    a(" -Ou | ",bcftools_bin," query -f '",query,"'")

  # Optional output file 
  if(!is.null(output_name))
    a(" -o ",output_name)

  if(!is.null(output_type))
    a(" -O",output_type[1])

  if(!is.null(extra_cmds))
    a(" ",extra_cmds)

  # Notes to self:
  # Why change regions contig names to match the VCF and not vice versa? I.e. why not "bcftools annotate --rename-chrs" before giving "-R"?
    # B/c "bcftools concat" must go first since only it can accept multiple files, and "-R" must be given to "concat" and can't be delayed
    # b/c after "concat" the file must be re-indexed, and "-R" requires an indexed file to work.
    # The file can't be indexed without writing an intermediate vcf/bcf file and running a separate "bcftools index" command.
    # Since this function's goal is to return just a single, pipable bcftools_cmd for seamless use, I can't do this.
  # Why -Ou before every pipe?
    # Speed. See bcftools docs for -O option and also https://www.biostars.org/p/336800/#336873.

  return(bcftools_cmd)
}
