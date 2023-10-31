#!/bin/bash
#
# This script converts 6 field BED files to GTF format


if [ x"$1" == x ]; then

        echo "please specify the name of the input BED file"

        exit 1

fi

if [ x"$2" == x ]; then

        echo "please specify the name of the output GTF file"

        exit 1

fi


awk '{OFS="\t";num++;print $1,"peak_"num, "exon", $2+1,$3,".", "+", ".", "gene_id peak_"num"; gene_name peak_"num";"}' $1 > $2
