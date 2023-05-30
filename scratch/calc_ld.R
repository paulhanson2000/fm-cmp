# TODO: Make so that user can supply reference file in any of vcf, gds, SeqArray gds, or PLINK formats.
# TODO: Add option to select by chrs and positions too. Maybe GRanges support too.
# TODO: Test if LDSTORE or PLINK are faster. If so, then just use that :V
  # Also ensure the different LD-calc methods give mostly identical results.
  # Or, do seqApply() + seqParallel() + C wizardry for max speed.

library(data.table)
library(SeqArray)
library(SNPRelate)

# Signed Pearson correlation.
# Returned LD matrix's rows/cols are sorted to match the order of the given rs_ids.
calc_ld <- function(ref_file, rs_ids, sample_ids, outfile_name=NULL, ref_alleles=NULL,
                    tmpfile_name="/tmp/tmp_snpgds_format_file.gds") {
  if(is.null(ref_alleles)) warning("calc_ld() warning: specifying your data's ref_alleles is highly recommended! If the reference alleles (i.e. non-effect alleles) in your data =/= those in the reference panel, the LD will be incorrect.")

  ref_gds <- seqOpen(ref_file, allow.duplicate=T)
  seqSetFilterAnnotID(ref_gds, rs_ids)
  seqSetFilter(ref_gds, sample.id=sample_ids)

  # Gotcha: SeqArray GDS and SNP GDS are in fact different formats, despite both using the .gds file extension.
  # snpgdsLDMat() must be used on an SNPGDSFileClass object, or will give faulty output.
  seqGDS2SNP(ref_gds, tmpfile_name, compress.geno="", compress.annotation="")
  tmpfile <- snpgdsOpen(tmpfile_name, readonly=F)
  snpgdsAlleleSwitch(tmpfile, toupper(ref_alleles))
  ld <- snpgdsLDMat(tmpfile, slide=0, method="corr", num.thread=getDTthreads())$LD

  # Sort LD matrix to match the order of the user-given rs_ids.
  snv_order <- match(rs_ids, seqGetData(ref_gds,"annotation/id"))
  ld <- ld[snv_order,snv_order]

  if(!is.null(outfile_name)) fwrite(ld, outfile_name, sep=' ', col.names=F)

  seqClose(ref_file)
  snpgdsClose(tmpfile)
  unlink(tmpfile_name)
  return(ld)
}
