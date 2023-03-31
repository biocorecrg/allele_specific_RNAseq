#!/bin/env bash

# Sum each allele counts for each file
# Author: Niccolo Arecco

if [ $# == 0 ]; then
    echo "Usage: $0 ARG1 "
    echo -e "\t- ARG1: <PATH/TO/output/cut_{1,2}/Allele_merged_Counts>"
    exit 1
fi

# Define input dir as $1
ALLELE_COUNTS_DIR="$1"
ALLELE_COUNTS_FILE=$( ls -1 ${ALLELE_COUNTS_DIR}/*.counts)

# For each file in the folder
for file in ${ALLELE_COUNTS_FILE} ; do
    ALLELE_COUNTS_FILENAME=$( basename $file .counts)
    ALLELE_COUNTS_FILENAME=$( echo ${ALLELE_COUNTS_FILENAME} | sed 's/_group//' )
    
    # Read the counts filename and the corresponding file
    # Skip the first line
    # Sum all the numbers in columns 2-5 without signif digits 
    # Output in a tab delimited fashion
   
    awk '{OFS = "\t"} NR>=2 { sA += $2 ; sB += $3 ; sU += $4 ; sAM += $5 } END \
        { print n1, sprintf("%.0f", sA), sprintf("%.0f", sB), sprintf("%.0f", sU), sprintf("%.0f", sAM) }' \
        n1=$ALLELE_COUNTS_FILENAME ${file}
done

