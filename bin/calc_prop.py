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
	usage = "usage: %prog -c <COMPOSITE COUNT TAB> -a <ALLELE COUNT TAB> -o <OUTPUT FILE>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-c', '--composite-input', help='Input composite count file', dest="starcount" )
	parser.add_option('-a', '--allele-input', help='Input allele count file', dest="allelecount" )
	parser.add_option('-s', '--strand', help='Input strandness (f/r/u)', dest="strand" )
	parser.add_option('-g', '--gtf', help='GTF file with annotation', dest="gtf" )
	parser.add_option('-o', '--output',help='ouput File', dest="wotus" )
	(opts,args) = parser.parse_args()
	if opts.starcount and opts.gtf and opts.allelecount and opts.wotus:pass
	if opts.strand == "f" or opts.strand == "r" or opts.strand == "u":pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	star_res = parseStar(opts.starcount, opts.strand)
	allele_res = parseAllele(opts.allelecount)
	annotations = parseGTF(opts.gtf)
	combineRes(star_res,allele_res,annotations,opts.wotus)


#AUXILIAR MODULES
def parseStar(file, strand):
	gene_comp_counts = collections.defaultdict(dict) 
	pp = pprint.PrettyPrinter(indent=4)
	fr = open(file ,'r')
	for line in fr:
		line = line.rstrip()
		fields = line.split("\t")
		if (line[:2]!="N_"):
			field = 0
			if (strand == "f"):
				field = 2
			if (strand == "r"):
				field = 3
			if (strand == "u"):
				field = 1
			geneid = fields[0]
			gene_count = fields[field]
			gene_comp_counts[geneid]	= gene_count
	return(gene_comp_counts)

def parseAllele(file):
	gene_all_counts = collections.defaultdict(dict) 
	pp = pprint.PrettyPrinter(indent=4)
	fr = open(file ,'r')
	for line in fr:
		line = line.rstrip()
		fields = line.split("\t")
		geneid = fields[0]
		allA_count = fields[1]
		allB_count = fields[2]
		gene_all_counts[geneid]["A"]	= allA_count
		gene_all_counts[geneid]["B"]	= allB_count
	return(gene_all_counts)

def parseGTF(file):
	gene_names = collections.defaultdict(dict) 
	pp = pprint.PrettyPrinter(indent=4)
	fr = open(file ,'r')
	for line in fr:
		line = line.rstrip()
		if (line[0]!="#"):
			fields = line.split("\t")
			if (fields[2] == "exon"):
				desc_dict = collections.defaultdict(dict)
				desc_fields = fields[8].rstrip(";").split(";")
				for desc_field in desc_fields:
					desc_vals = desc_field.lstrip().split(" ")
					desc_dict[desc_vals[0]] = desc_vals[1]
				gene_id = desc_dict["gene_id"].strip('"')
				if ("gene_name" in desc_dict):
					gene_name = desc_dict["gene_name"].strip('"')
				else:
					gene_name = gene_id
				gene_names[gene_id]["name"] = gene_name
				gene_names[gene_id]["chr"] = fields[0]
	return(gene_names)

def combineRes(count_res,allele_res,genames,wotus):
	f = open(wotus,'w')
	f.write("gene_id" + "\tgene_name\tchromosome\t" + "A/(A+B)" + "\t" + "B/(A+B)" + "\t" + "A/(A+B)*Comp" + "\t" + "B/(A+B)*Comp" + "\n")
	for geneid in count_res:
		comp_count = int(count_res[geneid])
		gene_name = genames[geneid]["name"]
		gene_chr = genames[geneid]["chr"]
		all_A = int(allele_res[geneid]["A"])
		all_B = int(allele_res[geneid]["B"])
		if (all_A == 0 and all_B == 0):
			propAB = 0
			propBA = 0
			calcAB = 0
			calcBA = 0
		else:
			propAB =  all_A/(all_A+all_B)
			propBA =  all_B/(all_A+all_B)
			calcAB = round((propAB*comp_count), 4)
			calcBA = round((propBA*comp_count), 4)
		f.write(geneid + "\t" + gene_name + "\t" + gene_chr + "\t" + str(propAB) + "\t" + str(propBA) + "\t" + str(calcAB) + "\t" + str(calcBA) + "\n")
	f.close()
	

#Calling
opts = options_arg()
if opts.starcount and opts.gtf and opts.allelecount and opts.wotus:
	__main__()

