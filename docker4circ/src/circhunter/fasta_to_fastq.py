#!/usr/bin/env python

import os, sys

"""
FASTA to FASTQ converter
"""
# Opening input file
input_data = open(sys.argv[1]).readlines()


new_object = ["",""]
circ_dict = {}

"""
FUNCTION: converter
OBJECTIVE: convert FASTA to FASTQ

x: refers to data
n: refers to new_object
"""


def converter(x, n):
	
	# Defining key of dictionary entry
	if x[0] == ">":
		n[1] = str()
		n[0] = ''.join(x[1:-1])
		return n;
	
	# Defining argument of dictionary entry	
	else:
		circ_partial = ''.join(x).strip()
		n[1] = n[1] + circ_partial
		return n;
		

# Main cycle
for line in input_data:
	data = list(line)
	converter(data, new_object)
	circ_dict[new_object[0]] = new_object[1]
	
for key in circ_dict:
	print("@%s" % key)
	print(circ_dict[key])
	print("+")
	
	score_length = len(circ_dict[key])
	print("I" * score_length)
