#!/bin/bash

# Special hardcoded case just for me and friends using Compute Canada, the data is already here. Just link to it instead of downloading.
if [ -d /lustre03/project/6074892/datasets/ ]; then
  ln -s /lustre03/project/6074892/datasets/ data
  echo "Success! The data can be found in the \"data/\" symbolic link that has been added to this directory."
fi

# Mahajan 2022 summary stats and fine-mapping results
mkdir -p data/DIAMANTE2022/sumstat/
cd       data/DIAMANTE2022/sumstat/
if ! [ -f DIAMANTE-SAS.sumstat.txt ]; then 
  curl -O https://personal.broadinstitute.org/ryank/DIAMANTE.sumstats.zip # Alternative: http://diagram-consortium.org
  unzip DIAMANTE.sumstats.zip
  parallel -t -j4 gzip -d ::: DIAMANTE-EAS.sumstat.txt.gz \
                              DIAMANTE-EUR.sumstat.txt.gz \
                              DIAMANTE-SAS.sumstat.txt.gz \
                              DIAMANTE-TA.sumstat.txt.gz
  rm DIAMANTE.sumstats.zip
fi
cd -

cd data/DIAMANTE2022/
if ! [ -d fine_mapping/ ]; then
  curl -O https://personal.broadinstitute.org/ryank/DIAMANTE.fine_mapping.zip
  unzip DIAMANTE.fine_mapping.zip
  rm    DIAMENTE.fine_mapping.zip
  mv fine_mapping_upload/ fine_mapping/
fi
cd -

# Reference Panels
# Only the reference panels' sample info is downloaded.
# Most data is downloaded later, in the R script, so that we only download what variants' data we actually need.
## 1000 Genomes reference panel sample info
mkdir -p data/ref/1kg/sample_info/
cd       data/ref/1kg/sample_info/
if ! [ -f integrated_call_samples_v3.20130502.ALL.panel ]; then
  curl -O http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
fi
cd -

## gnomAD 1000 Genomes + Human Genome Diversity Project reference panel sample info
mkdir -p data/ref/gnomad_1kg+hgdp/sample_info/
cd       data/ref/gnomad_1kg+hgdp/sample_info/
if ! [ -f gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv ]; then
  curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz
  gunzip -c gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz > gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv 
  rm gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz
fi
cd -

# Annotations
## Pancreatic islet enhancers enriched in T2D (Pasquali et al. 2014)
mkdir -p data/anno/pasquali/
cd       data/anno/pasquali/
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
mkdir -p data/anno/miguel-escalada/
cd       data/anno/miguel-escalada/
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
                               #https://www.crg.eu/sites/default/files/crg/islet_chromhmm.bed.zip
                               #https://www.crg.eu/sites/default/files/crg/pi_merged_washu_text.txt
                               #https://www.crg.eu/sites/default/files/crg/islet_regulome_simplified.bed.zip
fi
cd -

## PAINTOR-curated annotation library
### TODO: There is no way to automate the download AFAIK, the Dropbox uses temporary download links. So... eventually need derive the annotations from their original sources to make this automatic
### The below does not work.
#mkdir -p data/anno/paintor_anno_lib/
#cd       data/anno/paintor_anno_lib/
#if ! [ -f <><><><> ]; then
#  curl -O 'https://public.boxcloud.com/d/1/b1!2VWZxP2FVqFHWephuFlmRN9z1ijQGdY2iSIdZLM1UZpYAT5jP-_uOXzltWJxUcKBFuySrH9dFP-Lky-ViDN3hJVZsStGxyBDnfnLTRVGlFq7EnkjYGv4oDkt8dy1AtxCXmSVOFC51UMsdOwOWqH_EV36sW1B06NhM-DwYnsS59ASY_8scltmoYquasWzH7EuCLqggN74q2jVhpWkjv8XGshWSh-LaX8FgWRykYxdtrNZ_ssebZhB-2rOZsLUwKRVpc7vHT2jRPRITxUZPQxp7o1JCAHcRDpQ4bEASmFV-LWNkqUlngDGk4QkRAce-F0ezz-vCU4UlZDJX7K20KjMbDOkELGRf_58oH2wdPXEbWKX3b2DH3Hut5_Xp-6RejL080RoRODlPIS06sPT7vondxBf1LAOTP7PMDebUJX0fJBAcePBQV_qSWEvVle7qx7Ek5i_oCYX2Mxv_CsLYOFdC_kI6VZjOCSsIF0dbWC8lJZYakjIhXKr06zMkcypQGzqnuGfRuJ_WlrHWfnAFu5rd91BjUBQXNumi5y3uO5uh3vKho8uRLwZOYlw213wef0h0z4IwdAEuhbXdsfHMLrWaJQxSr94SzuX48OiYC2kan6mg0KPR91rOXiuKhZkh6l5PAwbQP-T9f66_8V2-IrNBa0udpYyJm5X5m4AeGF313M5z4waYuwzCSCMmfgup1t9Dfp8hkChUxLI54itHuDKRUSwRAWQEkj_v8KRZd_rr-XzedRTBPGkFRVaJ_jk4YzovY7vXQVG6umc_5V6JobEECMgVCCLJYCl8PHff7KIJazPcy2abQJG2_af3JNJPvsmWQ-y3c_gT7B36LN3a98jgLISMNUdPF6y4tuFyrp8DORtBJzCE9lbTf1ebsKsc_SEGP5I5kKK8ClEDJQ41YbxnA1zJTmJKvg4N3s6km8KFVgE6mrcMY0upggl2e_uI37crpRHQeRNs2BbvMnjWd2NVc03WXF_HejDrEQtxUhTVvrPLAxOiz3rxdtskJBqHETlKYri1bxQX5_SaI1AVhusmyKzlqr2hKEYhTrAefIJMyeNnAGAF1Q7lWNhQ-w67E7Hu88zfs9MKtzGwud0eYALlUEA4W7xrhUrsAMViM476ZwfBay5gFeKEB-w5OKbWlXk0cYmuLVKMosO9pfxp5lh5L481O3-w8Z2-5o-RS-hJn3D5Ri1Qt4LOmkZ5huiDTv2or4hD60SH-khihQbR8PYLIvNVsPQOXXVr82v1JO2sXfFgdeC3hSj2mvrCKPfxWem7Cb33_gR8QSMso-K6d0./download'
#fi
