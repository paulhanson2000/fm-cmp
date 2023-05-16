#!/bin/bash

# module purge # Good practice to do this, but comment out for a speedup
module load python/3.11.2
module load scipy-stack
module load plink/1.9b_6.21-x86_64 # `module load plink/2.00a3.6` (plink 2.0) doesn't seem to work

sumstats_eas="${SCRATCH}/DIAMANTE-EAS.sumstat-for_susiex.txt"
sumstats_eur="${SCRATCH}/DIAMANTE-EUR.sumstat-for_susiex.txt"
sumstats_sas="${SCRATCH}/DIAMANTE-SAS.sumstat-for_susiex.txt"
cat "../data/DIAMANTE2022/sumstat/DIAMANTE-EAS.sumstat.txt" | sed -e '1 s/Fixed-effects_beta/BETA/' | grep -Eiv "NA|NaN" > "${sumstats_eas}"
cat "../data/DIAMANTE2022/sumstat/DIAMANTE-EUR.sumstat.txt" | sed -e '1 s/Fixed-effects_beta/BETA/' | grep -Eiv "NA|NaN" > "${sumstats_eur}"
cat "../data/DIAMANTE2022/sumstat/DIAMANTE-SAS.sumstat.txt" | sed -e '1 s/Fixed-effects_beta/BETA/' | grep -Eiv "NA|NaN" > "${sumstats_sas}"

ld_eas="${SCRATCH}/susiex_ld_eas"
ld_eur="${SCRATCH}/susiex_ld_eur"
ld_sas="${SCRATCH}/susiex_ld_sas"

ref_eas=../data/reference_panels/1000_genomes/eas/g1000_eas
ref_eur=../data/reference_panels/1000_genomes/eur/g1000_eur
ref_sas=../data/reference_panels/1000_genomes/sas/g1000_sas

out_dir=../scratch/susiex_out

python3 ../third_party/SuSiEx/SuSiEx.py \
  --sst_file=${sumstats_eas},/scratch/phanso1/DIAMANTE-EUR.sumstat-for_susiex.txt,/scratch/phanso1/DIAMANTE-SAS.sumstat-for_susiex.txt \
  --n_gwas=39,32,15 \
  --ld_file=${ld_eas},${ld_eur},${ld_sas} `# If pre-existing file(s), computation will be skipped ` \
  --ref_file=${ref_eas},${ref_eur},${ref_sas} \
  --plink=plink \
  --out_dir=./susiex_out/ \
  --out_name=test \
  --chr=2 `#TODO` \
  --bp=27230940,28230940 `#TODO` \
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
  --keep-ambig=True
  # defaults: --max_iter=100, --pval_thresh=1e-5, --tol=1e-4, --n_sig=5, --level=95%, --min_purity=0.5, --full_out=NULL


rm ${sumstats_eas} ${sumstats_eur} ${sumstats_sas}
#rm ${SCRATCH}/susiex_ld*
