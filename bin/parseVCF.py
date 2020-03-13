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
import gzip 
from collections import Counter


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -i <input vcf file> -o <output vcf file>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-i', '--input', help='Input vcf file', dest="input" )
	parser.add_option('-1', '--species1', help='Input species id 1', dest="species1" )
	parser.add_option('-2', '--species2', help='Input species id 2', dest="species2" )
	parser.add_option('-o', '--output',help='ouput vcf File', dest="wotus" )
	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	parsefile(opts.input, opts.wotus, opts.species1, opts.species2)
	#longgenes = getlongest(myfastas)
	#write_seqs(myfastas, longgenes, opts.wotus)


#AUXILIAR MODULES
def parsefile(file, ofile, sp1, sp2):
	f = open(ofile,'w')
	pp = pprint.PrettyPrinter(indent=4)
	fieldnames = {}
	# prepare files files
	infile = open(file, 'r')
	if (file.endswith('.gz')):
		infile = gzip.open(file, 'rt')
	for line in infile:
		line = line.rstrip()
		if (line[:2]=="##"):
			f.write(line + "\n")
		elif(line[0]=="#"):
			fields = line.split('\t')
			ind1 = fields.index(sp1)
			ind2 = fields.index(sp2)
			del(fields[10:])
			fields[9] = sp1 + "|" + sp2
			line = "\t".join(fields)
			f.write(line + "\n")
		else:
			fields = line.split('\t')
			genA = fields[ind1].split(":", 1)[0]
			genB = fields[ind2].split(":", 1)[0]
			# Exclude cases in which on one of two alleles there are no genotype calls
			if (genA == './.' or genB == './.'):
				continue
			# Exclude cases in which on both calls are identical to the reference
			elif(genA == '0/0' and genB == '0/0'):
				continue
			else:
				newGenA = parseGen(genA) 
				newGenB = parseGen(genB) 
				if (newGenA >= 0 and newGenB >= 0 and fields[6] == "PASS"):
					del(fields[10:])
					fields[9] = str(newGenA) + "|" + str(newGenB)
					line = "\t".join(fields)
					f.write(line + "\n")
	infile.close()
	f.close()
	return


def parseGen(genotype):
	newGen = -1
	gen = int(genotype.split("/")[0])
	if (gen == 0):
		newGen = 0
	elif (gen < 3):
		newGen = 1
	return newGen
	
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

