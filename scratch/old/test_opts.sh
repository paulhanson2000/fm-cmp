#!/bin/bash

# getopt is not ideal. If I were reading this code myself I would have no idea what it means. So tossing this script.

ALPHA=unset
BETA=unset
CHARLIE=unset
DELTA=unset

usage() {
  echo -e "Usage: $0 [-i|--interactive] --:sumstat_file 
    -i, --interactive \t\t <>
    --sumstat_file \t\t (required)"
  exit 2
}

# 
PARSED_ARGS=$(getopt -n $0 --long interactive,\
sumstat_file:,\
ref_file:,\
sumstat_n,\
rsids,\
chr,\
pos_min,\
pos_max,\
rsid_col:,chr_col:,pos_col:,b_col:,se_col:,ref_allele_col: -- "$@")
ALL_ARGS_VALID=$?
if [ "$ALL_ARGS_VALID" != "0" ]; then usage; exit 1; fi
echo "PARSED_ARGS is $PARSED_ARGS"
eval set -- "$PARSED_ARGS" # Apparently eval bad, fix later. https://stackoverflow.com/questions/35235707/bash-how-to-avoid-command-eval-set-evaluating-variables

while :
do
  case "$1" in
    --interactive) interactive=1      ; shift  ;;
    --sumstat_file)  sumstat_file="$2"     ; shift 2;;
    --ref_file) ref_file="$2" ; shift 2;;

    --) shift; break ;; # "--)" means the end of the arguments; drop this, and break out of the while loop.
    *) echo "Unexpected option: $1 - this should not happen."; usage;; # If invalid options were passed, then getopt should have reported an error, which we checked as VALID_ARGUMENTS when getopt was called...
  esac
done

echo "ALPHA   : $sumstat_file"
echo "BETA    : $ref_file "
echo "CHARLIE : $interactive"
echo "Parameters remaining are: $@"
