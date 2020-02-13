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
	parser.add_option('-o', '--output',help='ouput prefix File', dest="wotus" )
	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	myreads = parsefile(opts.input, opts.wotus)


#AUXILIAR MODULES
def parsefile(file, oprefix):
	import pysam
	pp = pprint.PrettyPrinter(indent=4)
    # prepare out files
	samfile = pysam.AlignmentFile(file, "rb")
	ofiles = [oprefix + "_reference.bam", oprefix + "_alternative.bam", oprefix + "_ambiguous.bam"]
	pyoffiles = []
	for ofile in ofiles:
		pyoffiles.append(pysam.AlignmentFile(ofile, "wb", template=samfile))
	firstStatus = 0
	firstPair = ""
	readnum = 0
	status = 0
	for read in samfile.fetch(until_eof=True):
		# If the WASP filtering is passed:
		if(read.has_tag("vW")):
			if (read.get_tag("vW") == 1):
				readnum = readnum +1
				# Status 0 is reference, 1 is alternative, 2 is ambiguous
				# no need for more checks if the first pair is ambiguous
				if firstStatus != 2:
					alleles = read.get_tag("vA").tolist()
					alleles_count = Counter(alleles)
					# If reference (no variant) is found
					if 3 in alleles_count.keys():
						status = 0
					# If both alleles are found (i.e. it is ambiguous)
						if 1 in alleles_count.keys():
							status = 2
					# If only allele two is found
					elif 1 in alleles_count.keys():
						status = 1
					# If pairs are assigned to different alleles they are ambiguous
					if status != firstStatus and readnum == 2:
						status = 2
#				print(read) 
#				print(status)
			# reset after second pair info and write the results
			if (readnum == 2):
				pyoffiles[status].write(firstPair)
				pyoffiles[status].write(read)
				firstPair = ""
				firstStatus = 0
				readnum = 0
				status = 0
			else:
				firstPair = read
				firstStatus = status
	return



	

#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()


