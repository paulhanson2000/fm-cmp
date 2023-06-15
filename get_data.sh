#!/bin/bash

# Special case for friends using Compute Canada, the data is already here. Just link to it instead of downloading.
if [ -d /lustre03/project/6074892/datasets/ ]; then
  ln -s /lustre03/project/6074892/datasets/ data
  echo "Success! The data can be found in the \"data/\" symbolic link that has been added to this directory."
  exit 0
fi

# DIAMANTE 2022 summary stats and fine-mapping results
mkdir -p ./data/DIAMANTE2022/sumstat/
cd       ./data/DIAMANTE2022/sumstat/
wget https://personal.broadinstitute.org/ryank/DIAMANTE.sumstats.zip # Alternative: http://diagram-consortium.org
unzip DIAMANTE.sumstats.zip
gzip -d DIAMANTE-EAS.sumstat.txt.gz DIAMANTE-EUR.sumstat.txt.gz DIAMANTE-SAS.sumstat.txt.gz

cd -

cd ./data/DIAMANTE2022/
wget https://personal.broadinstitute.org/ryank/DIAMANTE.fine_mapping.zip
unzip DIAMANTE.fine_mapping.zip
rm    DIAMENTE.fine_mapping.zip
mv fine_mapping_upload/ fine_mapping/

cd -

# 1000 Genomes reference
## Sample information (gender, ancestry)
mkdir -p ./data/ref/1kg/sample_info/
cd       ./data/ref/1kg/sample_info/
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

cd -

## Genetic information
## TODO: Make the PLINK format myself by converting from the GDS format, instead of separate download?
mkdir -p ./data/ref/1kg/gds_format/
cd       ./data/ref/1kg/gds_format/
wget https://gds-stat.s3.amazonaws.com/download/1000g/2013/1KG_ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds

cd -

mkdir -p ./data/ref/1kg/plink_format/
cd       ./data/ref/1kg/plink_format/
parallel -t -j 3 wget  ::: https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip \
                           https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip \
                           https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_sas.zip
parallel -t -j 3 unzip ::: g1000_eas.zip \
                           g1000_eur.zip \
                           g1000_sas.zip
mkdir eas eur sas
mv g1000_eas* eas/
mv g1000_eur* eur/
mv g1000_sas* sas/

cd -
