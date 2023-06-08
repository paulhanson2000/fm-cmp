#!/bin/bash

# module purge # Good practice to do this, but comment out for a speedup
module load python/3.11.2
module load scipy-stack
module load plink/1.9b_6.21-x86_64 # `module load plink/2.00a3.6` (plink 2.0) doesn't seem to work

sumstats_eas=${SCRATCH}/DIAMANTE-EAS.sumstat-for_susiex.txt
sumstats_eur=${SCRATCH}/DIAMANTE-EUR.sumstat-for_susiex.txt
sumstats_sas=${SCRATCH}/DIAMANTE-SAS.sumstat-for_susiex.txt
cat "../../data/DIAMANTE2022/sumstat/DIAMANTE-EAS.sumstat.txt" | sed -e '1 s/Fixed-effects_beta/BETA/' | grep -Eiv "NA|NaN" | awk '$5=toupper($5)' | awk '$6=toupper($6)' > ${sumstats_eas}
cat "../../data/DIAMANTE2022/sumstat/DIAMANTE-EUR.sumstat.txt" | sed -e '1 s/Fixed-effects_beta/BETA/' | grep -Eiv "NA|NaN" | awk '$5=toupper($5)' | awk '$6=toupper($6)' > ${sumstats_eur}
cat "../../data/DIAMANTE2022/sumstat/DIAMANTE-SAS.sumstat.txt" | sed -e '1 s/Fixed-effects_beta/BETA/' | grep -Eiv "NA|NaN" | awk '$5=toupper($5)' | awk '$6=toupper($6)' > ${sumstats_sas}

ld_eas=${SCRATCH}/susiex_ld_eas
ld_eur=${SCRATCH}/susiex_ld_eur
ld_sas=${SCRATCH}/susiex_ld_sas

ref_eas=../../data/ref/1kg/plink_format/eas/g1000_eas
ref_eur=../../data/ref/1kg/plink_format/eur/g1000_eur
ref_sas=../../data/ref/1kg/plink_format/sas/g1000_sas

python3 ../../third_party/SuSiEx/SuSiEx.py \
  --sst_file=${sumstats_eas},${sumstats_eur},${sumstats_sas} \
  --n_gwas=56268,80154,16540 `#39,32,15` \
  --ld_file=${ld_eas},${ld_eur},${ld_sas} `# If pre-existing file(s), computation will be skipped ` \
  --ref_file=${ref_eas},${ref_eur},${ref_sas} \
  --plink=plink \
  --out_dir=out \
  --out_name=test \
  --chr=2 `#TODO` \
  --bp=27630940,27830940 \
  --chr_col=1,1,1 `# chromosome(b37) col` \
  --snp_col=4,4,4 `# rsID col` \
  --bp_col=2,2,2  `# position(b37) col` \
  --a1_col=5,5,5  `# effect_allele` \
  --a2_col=6,6,6  `# other_allele` \
  --eff_col=8,8,8 `# Fixed-effects_beta` \
  --se_col=9,9,9  `# Fixed-effects_SE` \
  --pval_col=10,10,10 `# Fixed-effects_p-value` \
  --mult-step=True \
  --maf=0.01 \
  --keep-ambig=True \
  --full_out=True `# output ALL snps used in fine-mapping in out.cs`
  # defaults: --max_iter=100, --pval_thresh=1e-5, --tol=1e-4, --n_sig=5, --level=95%, --min_purity=0.5


#rm ${sumstats_eas} ${sumstats_eur} ${sumstats_sas}
#rm ${SCRATCH}/susiex_ld*
