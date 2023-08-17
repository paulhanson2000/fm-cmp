Here is a rough explanation of the code structure to help understand it:

`init/` scripts just download data and build/install any needed software/packages.\
`run.qmd` runs the fine-mapping methods based on settings given in the config files.\
`compare\_results.qmd`'s job is to take the fine-mapping results and display them.

# `run.qmd`
## LiftOver
The given loci, sumstats, reference panel, and annotations must be lifted if their genomic coordinates are based on different reference genomes. It is better to avoid lifting the summary stats (and especially reference panel) if possible, as they are likely large files and lifting individual variants is more complicated than lifting regions. At minimum, we need:
+ Versions of the annotations for each of the sumstats' reference genome, so that annotation columns can be added to the sumstats later.
+ Loci lifted to the sumstats' reference genome(s) and reference panel's reference genome, because both the sumstats and panel need to be filtered by the loci regions.

If they are not already of the same build, we also require sumstats lifted to the reference panel's reference genome. This is done later in the script, after the summary stats and reference panel have been filtered a bit.\

## Loading data
Sumstats are loaded in and filtered down to only include the needed variants: those within the regions given in `loci.config` and with EAF>0.

Variants are *not* filtered by p-value, as some fine-mapping methods such as SuSiEx prefer this. If variants need to be filtered by p-value, it is done near the end in the final preparations for the individual methods.

Loading the reference panel is more complicated, as reference panels can be extremely huge (the gnomAD reference has files which are >100 GB *per chromosome*!). Therefore, combining bcftools and the `seqVCF2GDS` function from [SeqArray](https://bioconductor.org/packages/release/bioc/html/SeqArray.html), only the needed regions are piped over the internet directly into the highly storage-efficient SeqArray GDS format, without having to write any bulky intermediate files.

## rsIDs, "lifting" the sumstats
This step is only necessary if the sumstats and reference panel are of different reference genomes.\
Unfortunately, liftOver is not suitable for lifting individual variants ([source](https://genome.ucsc.edu/FAQ/FAQreleases.html#snpConversion)). Instead, it is recommended to lift variants over by their dbSNP rsID. Therefore, if the sumstats or reference panel does not have rsIDs, we search for them by their chromosome+position+alleles in dbSNP (dbSNP VCF file on the official ftp site [here](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/)).

In the end, the sumstats do not actually need to be lifted, we only need a set of shared IDs to be able to get information from the reference panel.\
There may be duplicated rsIDs or variants that could not be matched to an rsID. Barring a variant simply not having an assigned rsID yet, these duplicates/missing IDs are likely variants which have been merged/withdrawn in the newer versions of dbSNP. ((TODO) For now I just omit these, since it is only ~5% of them and it would be time-consuming to account for them, but it is possible: Ensembl's VEP can find old rsID synonyms. However, it takes an extremely long time to query. I am currently trying to find where I can download a list of all the archived dbSNP IDs for faster querying.)

## Filtering
We subset the sumstats to only the variants it shares with the reference panel. We cannot do anything with non-shared variants as they would be missing LD, so they are removed. MAC, MAF, and missingness filters are also applied (note the MAC, MAF, and missingness are calculated separately for each ancestry group.)

## Annotation
The sumstats are annotated by looking at the provided `.bed` files in `anno.config`. If a variant's position is within any of the regions provided in the `.bed` file, it is marked as positive for the annotation represented by that `.bed` file. The result is a bunch of columns of 0's and 1's being appended to the sumstats files.

`fgwas` is then run on these annotated sumstats to find which annotations are most significantly enriched in the most important variants. The result is, for each sumstat file, a list of the annotations which produced the highest-likelihood model, and prior probabilities for each variant based on this model.

(TODO) Currently working on how to incorporate these priors into some of the fine-mapping methods. As far as I know only susieR accepts custom priors as input.

## LD
For each combination of locus and ancestry group, an LD matrix is computed.

### Allele switching
There is no guarantee that the "effect"/"other" alleles in the sumstats are the same as the "reference"/"alternate" alleles in the reference panel. Therefore, the alleles must be switched. ((TODO) currently allele switching was applied to the reference panel just before calculating LD, but this must be changed to allow for pre-computed LD input. Instead, switching will be applied to the sumstats by flipping the betas/z-scores *after* ld calculation.)

## Running the fine-mapping methods
Each fine-mapping method has its own input requirements. The sumstats and LD are munged to fit these requirements, and the softwares are run.
