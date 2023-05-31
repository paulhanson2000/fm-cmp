library(data.table)
library(SeqArray)
library(SNPRelate)

# TODO: Make so that user can supply reference file in vcf, gds, SeqArray gds, or PLINK formats.
# TODO: add option to select by chrs and positions too. Maybe GRanges support too.

# TODO: make it so that the LD matrix is in the sam eorder as the rsids specified.

calc_ld <- function(ref_file, rs_ids, sample_ids, outfile_name=NULL, effect_alleles=NULL,
                    tmpfile_name="/tmp/tmp_snpgds_format_file.gds") {
  if(is.null(effect_alleles)) warning("calc_ld() warning: specifying effect_alleles is highly recommended! If the effect alleles in your data =/= the alternate alleles in the reference panel, the LD will be wrong.")
  ref_gds <- seqOpen(ref_file)

  seqSetFilterAnnotID(ref_gds, rs_ids)
  seqSetFilter(ref_gds, sample.id=sample_ids)
  seqGDS2SNP(ref_gds, tmpfile_name, compress.geno="", compress.annotation="")

  seqResetFilter(ref_gds)
  seqClose(reg_gds)

  tmpfile <- snpgdsOpen(tmpfile_name, readonly=F)
  snpgdsAlleleSwitch(tmpfile, toupper(effect_alleles))
  ld <- snpgdsLDMat(tmpfile, slide=0, method="corr", num.thread=getDTthreads())$LD

  # Sort LD matrix to match the order of the user-given rs_ids.
  snv_order <- match(rs_ids, seqGetData(ref_gds,"annotation/id"))
  ld <- ld[snv_order,snv_order]

  if(!is.null(outfile_name)) fwrite(ld, outfile_name, sep=' ', col.names=F)

  seqClose(ref_file)
  snpgdsClose(tmpfile)
  unlink(tmpfile_name)
  return(ld)

{seqSetFilter(ref_gds, variant.sel=1:10, sample.sel=1:2000)
rs_ids <- rev(seqGetData(ref_gds, "annotation/id"))
sample_ids <- seqGetData(ref_gds, "sample.id")
effect_alleles <- sapply(seqGetData(ref_gds, "allele"), function(x) strsplit(x,',')[[1]][2], USE.NAMES=F)
}
seqResetFilter(ref_gds)
calc_ld(ref_file = "../data/ref/1kg/gds_format/1KG_ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds",
        rs_ids = rs_ids,
        sample_ids = sample_ids,
        effect_alleles = effect_alleles)
  
ref_file = "../data/ref/1kg/gds_format/1KG_ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds"
tmpfile_name = "/tmp/tmp_snpgds_format_file.gds"

# Notes to self for peace of mind:
# Yes, read.gdsn(index.gdsn(tmpfile, "snp.rs.id")) is in the same order as seqGetData(ref_gsd,"annotation/id").
# No, GDS file does not care what order of rs_ids given in seqSetFilterAnnotID(), the rs_ids will be in w/e order they happen to be on disk.
# Yes, the SNPGDS tmpfile also retains w/e order the rs_ids happened to be in the SeqArray GDS, filtering never changes the order.
# snpgdsLDMat()$snp.id will be whatever the order is on disk too, no matter filters, no matter if you give snpgdsLDMAT() ~itself~ the snp.id arg neither. It doesn't care.
#
