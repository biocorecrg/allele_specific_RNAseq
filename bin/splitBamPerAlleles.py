#!/usr/bin/env python 
__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse
import collections
import pprint
import copy 
import subprocess
import textwrap 
from collections import Counter


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -i <input bam file> -o <OUTPUT PREFIX>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-i', '--input', help='Input bam file', dest="input" )
	parser.add_option('-m', '--minvar', help='Minimum number of variants', dest="minvar", default=1, type = int )
	parser.add_option('-o', '--output', help='Ouput prefix File', dest="wotus" )
	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	myfastas = parsefile(opts.input, opts.wotus, opts.minvar)


#AUXILIAR MODULES
def parsefile(file, oprefix, minvar):
	import pysam
	pp = pprint.PrettyPrinter(indent=4)
    # prepare out files
	samfile = pysam.AlignmentFile(file, "rb")
	ofiles = [oprefix + "_unknown.bam", oprefix + "_alleleA.bam", oprefix + "_alleleB.bam", oprefix + "_ambiguous.bam", oprefix + "_cutoff.bam"]
	pyoffiles = []
	for ofile in ofiles:
		pyoffiles.append(pysam.AlignmentFile(ofile, "wb", template=samfile))
	readnum = 0
	status = 0
	for read in samfile.fetch(until_eof=True):
		# If the WASP filtering is passed:
		if(read.has_tag("vW")):
			if (read.get_tag("vW") == 1):
				#readnum = readnum +1
				# the vA TAG gives a string with three possible values:
				# 1, 2, or three depending on how many SNPs and which genotype 
				# they match, are detected in each pair (so 1,1 means 2 SNPs 
				# matching genotype alle 1 considering both pairs) 
				alleles = read.get_tag("vA").tolist()
				alleles_count = Counter(alleles)
				status = 0;
				# If allele one is found
				if 1 in alleles_count.keys():
					status = 1
					# If both alleles are found (it is ambiguous)
					if 2 in alleles_count.keys():
						status = 3
				# If only allele two is found
				elif 2 in alleles_count.keys():
					status = 2
				#check the threshold
				if(status == 1):
					if(alleles_count[1] < minvar):
						status = 4
				if(status == 2):
					if(alleles_count[2] < minvar):
						status = 4
						
				pyoffiles[status].write(read)
	return

def count_allels(firstpair, secondpair):
	filtered_counts = dict()

	alleles_f = firstpair.get_tag("vA").tolist()
	alleles_s = secondpair.get_tag("vA").tolist()
	alleles_fc = Counter(alleles_f)
	alleles_sc = Counter(alleles_s)
	#for all_name in set(alleles_fc) | set(alleles_sc):
	#	filtered_counts[all_name]  = alleles_fc.get(all_name, 0) + alleles_sc.get(all_name, 0)
	#print(filtered_counts)
	
#	alleles_counts = Counter(alleles)
#	for variant_type in alleles_counts:
#		alleles_count = alleles_counts[variant_type]
#		if (alleles_count >= minvar):
#			filtered_counts[variant_type] = alleles_count
	return filtered_counts

	

	

#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()

