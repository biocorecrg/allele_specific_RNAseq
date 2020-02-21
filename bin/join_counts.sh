#!/bin/env bash
#
# This script will join differen read counts

if [ x"$1" == x ]; then

        echo "please specify file id"

        exit 1

fi

echo -e "gene_id\tgenoA\tgenoB\tref\tamb" > ${1}_group.counts
join ${1}_alleleA.counts ${1}_alleleB.counts | join - ${1}_ref.counts | join - ${1}_ambiguous.counts | tr ' ' '\t' >> ${1}_group.counts 

echo -e "#fname\tgenoA\tgenoB\tref\tamb" "" > ${1}_group.stats

awk -v fname=$1 '{if ($1!="_" || $1!="gene_id") {alla+=$2; allb+=$3; refs+=$4; ambs+=$5; }} END {print fname"\t"alla"\t"allb"\t"refs"\t"ambs}' ${1}_group.counts >> ${1}_group.stats

