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
	myfastas = parsefile(opts.input, opts.wotus)
	#longgenes = getlongest(myfastas)
	#write_seqs(myfastas, longgenes, opts.wotus)


#AUXILIAR MODULES
def parsefile(file, oprefix):
	import pysam
	pp = pprint.PrettyPrinter(indent=4)
    # prepare out files
	samfile = pysam.AlignmentFile(file, "rb")
	ofiles = [oprefix + "_undet.bam", oprefix + "_alleleA.bam", oprefix + "_alleleB.bam", oprefix + "_ambiguous.bam"]
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
				# No need for more checks if the first pair is ambiguous
				if firstStatus != 3:
					alleles = read.get_tag("vA").tolist()
					alleles_count = Counter(alleles)
					# If allele one is found
					if 1 in alleles_count.keys():
						status = 1
					# If both alleles are found (it is ambiguous)
						if 2 in alleles_count.keys():
							status = 3
					# If only allele two is found
					elif 2 in alleles_count.keys():
						status = 2
					# If pairs are assigned to different alleles they are ambiguous
					if status + firstStatus == 3:
						status = 3
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


def getlongest(seqs):
	longest_genes = collections.defaultdict(dict) 
	for geneid in seqs:
		gsize = 0
		mpepid = ""
		for pepid in seqs[geneid]:
			psize = len(seqs[geneid][pepid])
			if (psize > gsize):
				gsize = psize
				mpepid = pepid
		longest_genes[geneid] = mpepid
	return longest_genes
	
def write_seqs(seqs, long, ofile):
	f = open(ofile + "_longest.fa",'w')
	for geneid in long:
		protid = long[geneid]
		longseq = seqs[geneid][protid]
		fastseq =  ">" + protid + "|" + geneid + "\n" + textwrap.fill(longseq, 60) + "\n"
		f.write(fastseq)
	f.close()
	#print(longseq)
	

#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()


