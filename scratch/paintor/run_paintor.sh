#!/bin/bash
./PAINTOR -input input.files -in paintor_in/ -out paintor_out/ \
  -Zhead z_eas,z_eur,z_sas \
  -LDname ld_corr_eas,ld_corr_eur,ld_corr_sas \
  -annotations dummy_annot \
  -enumerate 1
