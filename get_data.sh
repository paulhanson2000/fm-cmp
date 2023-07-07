#!/bin/bash

# Special hardcoded case just for me and friends using Compute Canada, the data is already here. Just link to it instead of downloading.
if [ -d /lustre03/project/6074892/datasets/ ]; then
  ln -s /lustre03/project/6074892/datasets/ data
  echo "Success! The data can be found in the \"data/\" symbolic link that has been added to this directory."
fi

# Mahajan 2022 summary stats and fine-mapping results
mkdir -p ./data/DIAMANTE2022/sumstat/
cd       ./data/DIAMANTE2022/sumstat/
if ! [ -f DIAMANTE-EAS.sumstat.txt ]; then 
  wget https://personal.broadinstitute.org/ryank/DIAMANTE.sumstats.zip # Alternative: http://diagram-consortium.org
  unzip DIAMANTE.sumstats.zip
  parallel -t -j3 gzip -d ::: DIAMANTE-EAS.sumstat.txt.gz \
                              DIAMANTE-EUR.sumstat.txt.gz \
                              DIAMANTE-SAS.sumstat.txt.gz
fi
cd -

cd ./data/DIAMANTE2022/
if ! [ -d fine_mapping/ ]; then
  wget https://personal.broadinstitute.org/ryank/DIAMANTE.fine_mapping.zip
  unzip DIAMANTE.fine_mapping.zip
  rm    DIAMENTE.fine_mapping.zip
  mv fine_mapping_upload/ fine_mapping/
fi
cd -

# 1000 Genomes reference panel
## Sample information (gender, ancestry)
mkdir -p ./data/ref/1kg/sample_info/
cd       ./data/ref/1kg/sample_info/
if ! [ -f integrated_call_samples_v3.20130502.ALL.panel ]; then
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
fi
cd -

## Genetic information
## TODO: Make the PLINK format myself by converting from the GDS format, instead of separate download?
mkdir -p ./data/ref/1kg/gds_format/
cd       ./data/ref/1kg/gds_format/
if ! [ -f 1KG_ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds ]; then
  wget https://gds-stat.s3.amazonaws.com/download/1000g/2013/1KG_ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds
fi
cd -

mkdir -p ./data/ref/1kg/plink_format/
cd       ./data/ref/1kg/plink_format/
if ! [ -d sas/ ]; then
  parallel -t -j3 wget  ::: https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip \
                            https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip \
                            https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_sas.zip
  parallel -t -j3 unzip ::: g1000_eas.zip \
                            g1000_eur.zip \
                            g1000_sas.zip
  mkdir eas eur sas
  mv g1000_eas* eas/
  mv g1000_eur* eur/
  mv g1000_sas* sas/
fi
cd -

# TOPLD reference panel
## Variant annotation info, incluing ancestry-specific AFs.
mkdir -p ./data/ref/topmed/topld/anno/
cd       ./data/ref/topmed/topld/anno/
if ! [ -f SAS_topld_anno.csv ]; then
  parallel -t -j3 curl -O ::: http://topld.genetics.unc.edu/downloads/downloads/EAS/SNV/EAS_chr[1-22]_no_filter_0.2_1000000_info_annotation.csv.gz \
                              http://topld.genetics.unc.edu/downloads/downloads/EUR/SNV/EUR_chr[1-22]_no_filter_0.2_1000000_info_annotation.csv.gz \
                              http://topld.genetics.unc.edu/downloads/downloads/SAS/SNV/SAS_chr[1-22]_no_filter_0.2_1000000_info_annotation.csv.gz
  parallel -t -j3 gunzip -d ::: *
  # TODO: what about that "no preprocessing" philosophy for the data folder?
  head -1 EAS_chr1_* | tee -a EAS_topld_anno.csv EUR_topld_anno.csv SAS_topld_anno.csv # Header, then...
  awk 'FNR-1' EAS_chr* >> EAS_topld_anno.csv # ...combine files
  awk 'FNR-1' EUR_chr* >> EUR_topld_anno.csv
  awk 'FNR-1' SAS_chr* >> SAS_topld_anno.csv

  parallel -t -j
fi

## LD
mkdir -p ./data/ref/topmed/topld/ld/
cd       ./data/ref/topmed/topld/ld/
if ! [ -f SAS_chr22_no_filter_0.2_1000000_LD.csv.gz ]; then
  parallel -t -j3 curl -O ::: http://topld.genetics.unc.edu/downloads/downloads/EAS/SNV/EAS_chr[1-22]_no_filter_0.2_1000000_LD.csv.gz
                              http://topld.genetics.unc.edu/downloads/downloads/EUR/SNV/EUR_chr[1-22]_no_filter_0.2_1000000_LD.csv.gz
                              http://topld.genetics.unc.edu/downloads/downloads/SAS/SNV/SAS_chr[1-22]_no_filter_0.2_1000000_LD.csv.gz
fi
cd -


# Annotations
## Pancreatic islet enhancers enriched in T2D (Pasquali et al. 2014)
mkdir -p ./data/anno/pasquali/
cd       ./data/anno/pasquali/
if ! [ -f PDX1_HI_32_1e-10_peaks.bed ]; then
#wget -r -np -A bed ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/FOXA2_HI_32_1e-10_peaks.bed 
  parallel -t -j13 curl -O ::: ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/FOXA2_HI_101_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/H2AZ_HI_22_1e-5_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/H2AZ_HI_32_1e-5_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/H2AZ_HI_34_1e-5_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/MAFB_HI_81_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/MAFB_HI_87_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/NKX2_2_HI_87_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/NKX2_2_HI_88_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/NKX6_1_HI_102_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/NKX6_1_HI_118_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/PDX1_HI_45_1e-10_peaks.bed \
                               ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/919/E-MTAB-1919/Files/PDX1_HI_32_1e-10_peaks.bed
fi
cd -

## Pancreatic islet chromatin architecture stuff (Miguel-Escalada et al. 2019)
mkdir -p ./data/anno/miguel-escalada/
cd       ./data/anno/miguel-escalada/
if ! [ -f atac_consistent_peaks.bed ]; then
  # TODO: there were a couple other non-bed format files I skipped b/c I don't know how to use them
  parallel -t -j11 curl -O ::: https://www.crg.eu/sites/default/files/crg/smca1_consistent_peaks_q001.bed \
                               https://www.crg.eu/sites/default/files/crg/h3k4me3_consistent_peaks_q005.bed \
                               https://www.crg.eu/sites/default/files/crg/ctcf_consistent_peaks_q001.bed \
                               https://www.crg.eu/sites/default/files/crg/h3k27ac_consistent_peaks_q005.bed \
                               https://www.crg.eu/sites/default/files/crg/med1_consistent_peaks_q001.bed \
                               https://www.crg.eu/sites/default/files/crg/islet_enhancer_hubs.bed \
                               https://www.crg.eu/sites/default/files/crg/islet_pats.bed \
                               https://www.crg.eu/sites/default/files/crg/islet_super_enhancers.bed \
                               https://www.crg.eu/sites/default/files/crg/islet_tad-like_domains.bed \
                               https://www.crg.eu/sites/default/files/crg/atac_consistent_peaks.bed
fi
cd -
