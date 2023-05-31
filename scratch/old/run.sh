#!/bin/bash

# getopt is not ideal. If I were reading this code myself I would have no idea what it means. So tossing this script.

help() {
  echo "Usage: $0 "
}

options=$(getopt --longoptions sumstat_file:,ref_file:,out_file,sumstat_n <>)

while getopts ":h" option; do
  case $option in
    h) help; exit 1;;
    sumstat_file) sumstat_file=$OPTARG;;
    ref_file)

    out_file) out=$OPTARG;;

    sumstat_n) sumstat_n=$OPTARG;;
    rsids) rsids=$OPTARG;;
    chr) chr=$OPTARG;;
    pos_min) pos_min=$OPTARG;;
    pos_max) pos_max=$OPTARG;;

    rsid_col) rsid_col=$OPTARG;;
    chr_col) chr_col=$OPTARG;;
    pos_col) pos_col=$OPTARG;;
    b_col) b_col=$OPTARG;;
    se_col) se_col=$OPTARG;;
    ref_allele_col) ref_allele_col=$OPTARG;;

    *) echo "Invalid option or argument $OPTARG"; exit 1;;
  esac
done

echo "$sumstat_file"


# Parse options. Note that options may be followed by one colon to indicate 
# they have a required argument
if ! options=$(getopt -o abc: -l along,blong,clong: -- "$@")
then
    # Error, getopt will put out a message for us
    exit 1
fi

set -- $options

while [ $# -gt 0 ]
do
    # Consume next (1st) argument
    case $1 in
    -a|--along) 
      aflag="y" ;;
    -b|--blong) 
      bflag="y" ;;
    # Options with required arguments, an additional shift is required
    -c|--clong) 
      cargument="$2" ; shift;;
    (--) 
      shift; break;;
    (-*) 
      echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) 
      break;;
    esac
    # Fetch next argument as 1st
    shift
done
