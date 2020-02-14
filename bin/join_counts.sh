#!/bin/env bash
#
# This script will join differen read counts

if [ x"$1" == x ]; then

        echo "please specify file id"

        exit 1

fi

echo -e "gene_id\tref\talt\tamb" > ${1}_group.counts
join ${1}_reference.counts ${1}_alternative.counts| join - ${1}_ambiguous.counts | tr ' ' '\t' >> ${1}_group.counts 

echo -e "#fname\tgene_id\tref\talt\tamb" "" > ${1}_group.stats

awk -v fname=$1 '{if ($1!="_" || $1!="gene_id") {refs+=$2; alts+=$3; ambs+=$4}} END {print fname"\t"refs"\t"alts"\t"ambs}' ${1}_group.counts >> ${1}_group.stats

