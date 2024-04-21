risk1=$1
risk2=$(echo $1^2 | bc) # Risk squared
./hapgen2 -m data/example/ex.map -l data/example/ex.leg -h data/example/ex.haps -t data/example/ex.tags -o out -dl 1085679 1 1.5 2.25 2190692 0 $risk1 $risk2 -n 5000 5000 -no_haps_output #-output_snp_summary
# -no_haps-output b/c only need genotypes for association analysis.
