#!/bin/bash
risk1=$1
risk2=$(echo $1^2 | bc) # Risk squared
./hapgen2 \
  -m data/1kg/genetic_map_chr22_combined_b37.txt \
  -l data/1kg/chr22_mac_gt40.hap.gz \
  -h data/1kg/chr22_mac_gt40.leg.gz \
  -o out \
  -dl 9411485 1 $risk1 $risk2 \
  -n 10 10 \
  -no_haps_output #-output_snp_summary
# -no_haps_output b/c only need genotypes for association analysis.
