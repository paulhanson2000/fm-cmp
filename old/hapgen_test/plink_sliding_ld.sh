#!/bin/bash
# "module load StdEnv/2023 plink/2.00a5.8" is too old, need latest version from website.
#./plink2 --haps   data/1kg/1000GP_Phase3_chr22.hap.gz \
#         --legend data/1kg/1000GP_Phase3_chr22.legend.gz 22 \
#         --mac 40 minor \
#         --indep-pairphase 100kb 1 0.2
#./plink2 --haps   data/1kg/1000GP_Phase3_chr22.hap.gz \
#         --legend data/1kg/1000GP_Phase3_chr22.legend.gz 22 \
#         --mac 40 minor \
#         --r-phased \
#         --ld-window-kb 50 \
#         --ld-window-r2 0 # 0.2^2


# RM ME
../simgwas_test/plink2 \
  --haps   ../simgwas_test/data/chr22_mac_gt40.hap.gz \
  --legend ../simgwas_test/data/chr22_mac_gt40.leg.gz 22 \
  --mac 40 minor \
  --r-phased \
  --ld-window-kb 50 \
  --ld-window-r2 0 # 0.2^2
