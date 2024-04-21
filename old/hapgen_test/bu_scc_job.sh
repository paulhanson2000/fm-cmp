qsub -N 1_hapgen_test.sh -l h_rt=8:00:00 -pe omp 8 -l mem_per_core=8G -P t2d-intern
