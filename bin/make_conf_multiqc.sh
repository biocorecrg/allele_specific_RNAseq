#!/bin/env bash
#
# This script will 

if [ x"$1" == x ]; then
        echo "please specify a title file"
        exit 1
fi

if [ x"$2" == x ]; then
        echo "please specify a subtitle file"
        exit 1
fi

if [ x"$3" == x ]; then
        echo "please specify a PI"
        exit 1
fi

if [ x"$4" == x ]; then
        echo "please specify a user"
        exit 1
fi

if [ x"$5" == x ]; then
        echo "please specify a contact mail"
        exit 1
fi

if [ x"$6" == x ]; then
        echo "please specify a reference genome"
        exit 1
fi

if [ x"$7" == x ]; then
        echo "please specify a pre config file"
        exit 1
fi

cat > config.yaml << EOL
title: "$1"
subtitle: "$2"
intro_text: False

report_header_info:
- PI: $3
- User: $4
- Date: `date`  
- Contact E-mail: '$5'
- Application Type: 'Allele specific RNAseq'
- Reference Genome: $6
EOL

cat $7 >> config.yaml
