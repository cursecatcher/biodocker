#!/usr/bin/env python

import os,sys

"""
This script uses a biomart export to return all introns inferrable from
exons in transcripts.

NOTE
The export must be sorted in the following order:
ENSG, ENST, RANK
"""

file_dict = {} # Creating dictionary that will contain the input file

export_data = open(sys.argv[1]).readlines() # Opening input file
keygen = 0000000000 # Key generator for dictionary

"""
FUNCTION: find_intron
OBJECTIVE: to find an intron between two exons

x: refers to exon
y: refers to exon_next
"""

def find_intron(x,y):
	
	# checks to add intron: checks chr, ENST and rank
	if x[3] == y[3] and x[1] == y[1] and int(x[7]) < int(y[7]):
		
		# calculating intron values
		ENSG = x[0]
		ENST = x[1]
		ENSE = ("%si" % x[2])
		chromosome = x[3]
		strand = x[6]
		
		if strand == "1":
			start = str(int(x[5]) + 1)
			end = str(int(y[4]) - 1)
		elif strand == "-1":
			start = str(int(y[5]) + 1)
			end = str(int(x[4]) - 1)

		# Renaming ranks
		x_rank_raw = "00%s" % x[7]
		x_rank = x_rank_raw[-3:]
		intron_rank_raw = "00%si" % x[7]
		intron_rank = intron_rank_raw[-4:]
		
		# Defining intron data
		intron = [ENSG, ENST, ENSE,
		chromosome, start, end,
		strand, intron_rank, x[8],
		x[9], x[10], x[11]]
		
		# Redefining original exon, with renamed rank
		exon = [x[0], x[1], x[2],
		x[3], x[4], x[5],
		x[6], x_rank, x[8],
		x[9], x[10], x[11]]
		
		# Printing original exon and calculated intron
		print("%s" % ("\t".join(exon)))
		print("%s" % ("\t".join(intron)))
	
	else:
		# Renaming rank
		x_rank_raw = "00%s" % x[7]
		x_rank = x_rank_raw[-3:]
		
		# Redefining original exon
		exon = [x[0], x[1], x[2],
		x[3], x[4], x[5],
		x[6], x_rank, x[8],
		x[9], x[10], x[11]]
		
		# Printing exon
		print("%s" % ("\t".join(exon)))		
		


# Main cycle
# Reading input file
for line in export_data:
	data = line.strip().split("\t")
	file_dict[keygen] = data # Adding line to dictionary
	keygen += 1 # Increasing dictionary's key

"""
# Debug
for key in file_dict:
	print("%s\t%s" % (key, file_dict[key]))
"""

# Reading dictionary, two keys at a time
for key in file_dict:
	exon = file_dict[key] # Reading exon
	
	try:
		exon_next = file_dict[key + 1] # Reading next exon
		find_intron(exon, exon_next)
	except:
		# Renaming rank
		exon_rank_raw = "00%s" % exon[7]
		exon_rank = exon_rank_raw[-3:]
		
		# Redefining original exon
		fixed_exon = [exon[0], exon[1], exon[2],
		exon[3], exon[4], exon[5],
		exon[6], exon_rank, exon[8],
		exon[9], exon[10], exon[11]]
		
		# Printing exon
		print("%s" % ("\t".join(fixed_exon)))
		break
