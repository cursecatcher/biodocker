#!/usr/bin/env python

import os,sys

"""
This script calculates 70bp junction sequences of circRNAs without
considering the transcript (it utilizes circRNA coordinates only).
"""

circ_data = open(sys.argv[1]).readlines()


"""
FUNCTION: fake_exons
OBJECTIVE: obtain fake exons to use in hash of span +35/-35 bp
If the circRNA length is below the threshold (71), the function will
split the circRNA in half and obtain fake exons with the two halves.

x: assign to data
"""

def fake_exons(x):
	n = x[0].strip().split("_")
	
	# Defining circRNA features
	circ_name = x[0]
	chromosome = n[0]
	start = int(n[1])
	end = int(n[2])
	strand = x[1]
	length = end - start + 1

	# Checking threshold
	if length > 70:
		
		# Defining exons
		ex5_end = start + 34
		ex3_start = end - 34
		
		jlength = (ex5_end - start + 1) + (end - ex3_start + 1)
		
		output = [circ_name, chromosome, strand,
		str(start), str(ex5_end), str(ex3_start), str(end), str(jlength)]
		print("\t".join(output))
	

	else:
		check = end - start
		
		if check%2 == 0:
			center = (end - start) / 2
			
			ex5_start = start
			ex5_end = (start + center) - 1
			exon5 = [ex5_start, ex5_end]
			
			ex3_start = end - center
			ex3_end = end
			exon3 = [ex3_start, ex3_end]
			
			jlength = (ex5_end - start + 1) + (end - ex3_start + 1)
			
			output = [circ_name, chromosome, strand,
			str(ex5_start), str(ex5_end), str(ex3_start), str(ex3_end), str(jlength)]
			print("\t".join(output))
		
		else:
			center = (end - start - 1) / 2
			
			ex5_start = start
			ex5_end = start + center
			exon5 = [ex5_start, ex5_end]
			
			ex3_start = end - center
			ex3_end = end
			exon3 = [ex3_start, ex3_end]
			
			jlength = (ex5_end - start + 1) + (end - ex3_start + 1)
			
			output = [circ_name, chromosome, strand,
			str(ex5_start), str(ex5_end), str(ex3_start), str(ex3_end), str(jlength)]
			print("\t".join(output))




# Main cycle
for circ in circ_data:
	data = circ.strip().split("\t")
	fake_exons(data)
