# Config Options
## `data.config`
Tell the pipeline about your GWAS summary statistics files.
+ `filepath`
+ `n`: sample size.
+ `ancestry`: ancestry of the individuals that the summary stats were derived from. e.g. EAS, EUR, SAS.
+ `ref_genome`: reference genome the file's variant positions are based on (UCSC naming scheme, e.g. "hg19", "hg38")
+ `..._col`: specify which column number in the file mean what:
  + `chr_col`: chromosome
  + `pos_col`: position
  + `rs_id_col`: dbSNP rsID
  + `effect_allele_col`
  + `other_allele_col`
  + `eaf_col`: effect allele frequency
  + `b_col`: beta value
  + `se_col`: standard error
  + `z_col`: z-score
  + `p_col`: p-value
  + `n_col`: per-variant sample size

## `loci.config`
Configure the regions to be fine-mapped.
+ `locus`: locus name (TODO: '+' in locus names will cause fgwas bug)
+ `chr`, `pos_min`, `pos_max`: defines the locus's region
+ `ref_genome`: reference genome the locus's positions are based on (UCSC naming scheme, e.g. "hg19", "hg38")
Note: the same locus name can be used for multiple regions to fine map those regions together.

## `anno.config`
(Optional) Specify functional annotations in the form of [`.bed`-format files](http://www.genome.ucsc.edu/FAQ/FAQformat.html#format1).
+ `name`
+ `ref_genome`: reference genome the `.bed` file's positions are based on (UCSC naming scheme, e.g. "hg19", "hg38)
+ `bed_file`: path to the `.bed` file

## `misc.config`
Specify which reference panel to use, which fine-mapping methods to run, MAF/MAC thresholds, and more.
+ `scratch_dir`: directory to store intermediate outputs that are worth keeping so they don't have to be recomputed, such as LD.
+ `maf_threshold`: minimum minor allele frequency allowed in summary stats and reference panel. NULL to set no threshold.
+ `mac_threshold`: minumum minor allele count allowed in summary stats and reference panel. NULL to set no theshold.
+ `missingness_threshold`: maximum proportion of missing data variants in the reference panel are allowed to have. NULL to set no threshold.
+ `ref_panel`: reference panel used for LD and filtering. Many options:
  + 3 pre-implemented options are availabe, totally hands-free. All you need to do is specify one of the following strings:
    + `"gnomad_1kg+hgdp"`: (recommended) gnomAD combined 1000 Genomes 30x and Human Genome Diversity Project
    + `"1kg_30x"`: 1000 Genomes 30x coverage ([paper](https://doi.org/10.1016/j.cell.2022.08.004))
    + `"1kG"`: 1000 Genomes phase 3
  + VCF file(s) to be used as the reference panel. If you specify multiple files, make sure they have the same samples; for example, if you have files split by chromosome.
  + SeqArray GDS file
  + `"my_ld"`; will use pre-computed LD matrices specified in `config.ld`. Useful if you cannot have direct access to the individual-level data of the reference panel you wish to use.
+ ` gnomad_count_fin_as_eur`: true or false. The hands-free gnomAD reference panel has separate Finnush and non-Finnish European populations. Should these be combined into single EUR ancestry group?
+ `ref_panel_ref_genome`: if you are providing your own reference panel, specify which reference genome its positions are based on in UCSC naming scheme (e.g. "hg19", "hg39")
+ Sample info (if you are providing your own reference panel)
  + `sample_info_file`
  + `sample_id_field`: name of the sample ID column.
  + `sample_ancestry_field`: name of the sample ancestry column.

## `fm.config`
Specify which fine-mapping methods to run, what p-value thresholds to use for them, and if they should use built-in annotation-incorportation functionality if any.
+ `methods`: which fine-mapping methods to generate outputs for.
+ `pval_thresholds`: p-value threshold for each method. Some fine-mapping methods don't handle large amounts of variants well so this may be helpful. (TODO: list which methods this is/isn't recommended for.)
+ `use_anno_functionality`: (TODO) true or false for each method, whether to use built-in support for annotations or not. If true, will override other annotation methods such as fgwas.

## `ld.config` (TODO)
If you specified `misc_config$ref_panel="my_ld"`, use `config.ld` to specify your pre-computed LD files:
+ `filepath`: name of a space-delimited LD matrix file.
+ `locus`: name of the `config.loci` locus the LD matrix is for.
+ `ancestry`: ancestry group the LD was calculated from.
+ `colnames_file`: useful if your LD file is only numbers and corresponding variant IDs are stored in another file.

## `ucsc_chr_names_map.txt`
+ Specify which chromosome names (`from`) are synonymous with those under the UCSC naming scheme (`to`).
