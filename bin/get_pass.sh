#!/bin/env bash
#
# This script will cut the length of sequences within a fastq file

if [ x"$1" == x ]; then

        echo "please specify a snp vcf file (it can be either plain or gzipped)"

        exit 1

fi


if [ `echo $1 | grep "gz"` ]; then
	cat="zcat "
else 
	cat="cat "
fi


$cat $1 |awk '{if ($0~"#") {print} else {if ($7=="PASS") {print}} }' > pass.vcf
gzip pass.vcf
