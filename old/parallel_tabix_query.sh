#!/bin/bash
query_file=$1
n_thread=$2

#n_lines=$(wc -l < $query_file)
#echo $n_lines
#n_lines_per_split=$((n_lines/n_thread + 1))
#echo $n_lines_per_split
#split -l $n_lines_per_split $query_file

split -n l/$n_thread $query_file query_part_
parallel 'tabix GCF_000001405.40-rsid_reformatted.gz -R {} > query_result_{#}.txt' ::: query_part_*
cat query_result_* > query_result.txt
rm query_part_* query_result_*
