#!/usr/bin/env python

import os,sys

"""
LAST UPDATED VERSION
-> WRITING UNIVOCAL CLASSIFICATION

This is a script that allows to calculate the presumed linear sequences
that form between the exons near to the circRNA and a circRNA.

It takes as input:
- export of all the exons
- intersect of circRNAs with bedtools with those exons

Note: the script will consider only exons at the 5' and 3' of circRNAs
"""

# Defining dictionaries
exon_dic = {}
circ_dic = {}
sorted_circ_dic = {}
rank_dic = {}
sorted_rank_dic = {}
unaligned_exon_dic = {}
intertranscript_dic = {}

# Opening input files
exon_data = open(sys.argv[1]).readlines()
overlap_data = open(sys.argv[2]).readlines()


"""
FUNCTION: save_exons
OBJECTIVE: saves all exons from a biomart export into a dictionary

x: assign to data
d: assign to exon_dic
"""
def save_exons(x, d):	
	# Defining key
	exon_name = str("%s_%s" % (x[1], x[2]))
	
	if exon_name not in d:
		d[exon_name] = x

"""
FUNCTION: circ_transcripts
OBJECTIVE: to create a dictionary containing all the exons involved in
a circRNA

x: assign to data
d: assign to circ_dic
un: assign to unaligned_exon_dic
exons: assign to exon_dic
"""

def circ_transcripts (x, d, un, it, exons):
	# Defining key: composed by circ name and transcript (chr_start_end_ENST)
	n = x[3].strip().split("_")
	circ_name = str("%s_%s" % (x[8], n[0]))
	strand = x[4]
	
	if circ_name not in d:
		d[circ_name] = []
	
	# Defining ENSE
	e = x[3].strip().split("_")
	ENSE = e[1]
	
	# Defining the variables needed for exon checks
	ex_start = int(x[1])
	ex_end = int(x[2])
	circ_start = int(x[6])
	circ_end = int(x[7])
	
	# Searching for exon rank in exon_dic
	exon = x[3]
	rank = str(exons[exon][7])
	ENSG = exons[exon][0]
	ENST = exons[exon][1]
	
	# Creating tuple to append to dictionary
	new_entry = (rank, ENSE, strand, ENSG, ENST)
	
	# Checking exons and introns and appending ENSEs in the dictionary
	if ( (ex_start == circ_start and ex_end <= circ_end)
	or (ex_end == circ_end and ex_start >= circ_start) ):
		d[circ_name].append(new_entry)
	
	# extended check for introns overlapping start/end of circRNAs
	elif ((len(rank) == 4 and circ_start > ex_start and circ_start < ex_end)
	or (len(rank) == 4 and circ_end > ex_start and circ_end < ex_end)):
		d[circ_name].append(new_entry)
	
	# extended check for exons overlapping start/end of circRNAs to save
	# in the unaligned_exon_dic, for later putative_exon classification
	elif ((len(rank) == 3 and circ_start > ex_start and circ_start < ex_end)
	or (len(rank) == 3 and circ_end > ex_start and circ_end < ex_end)):
		
		if circ_name not in un:
			un[circ_name] = []
		
		un[circ_name].append(new_entry)
	
	# SECTION FOR INTERTRANSCRIPT circRNAs
	circ = x[8]
	
	if circ not in it:
		it[circ] = []
	
	# Appending first element in indiscriminate manner
	if len(it[circ]) == 0:
		if (ex_start == circ_start and ex_end <= circ_end):
			new = ("start", ENST)
			it[circ].append(new)
		
		elif (ex_end == circ_end and ex_start >= circ_start):
			new = ("end", ENST)
			it[circ].append(new)
	
	# Appending second element only if it is mapped in the opposite
	# side of the circRNA
	
	elif len(it[circ]) == 1:
		if (it[circ][0][0] != "start" and it[circ][0][1] != ENST
		and ex_start == circ_start and ex_end <= circ_end):
			new = ("start", ENST)
			it[circ].append(new)
		
		elif (it[circ][0][0] != "end" and it[circ][0][1] != ENST
		and ex_end == circ_end and ex_start >= circ_start):
			new = ("end", ENST)
			it[circ].append(new)
	
"""
FUNCTION: save_ranks
OBJECTIVE: create a dictionary filled with exon ranks, in order to later
calculate max and min rank of a transcript

x: refers to data
d: refers to rank_dic
"""

def save_ranks(x, d):
	
	ENST = x[1]
	rank = str(x[7])
	
	if ENST not in d:
		d[ENST] = []
	
	# Adding rank value to transcript in dictionary
	d[ENST].append((rank, x[2]))


"""
FUNCTION: seq5
OBJECTIVE: obtain the coordinates of the first exon of the circRNA and
the coordinates of the exon before it.

c: refers to circ
ex: refers to exon (obtained from sorted_circ_dic[circ][0])
rd: refers to sorted_rank_dic
ed: refers to exon_dic
l: refers to label
"""
"""
def seq5(c, ex, rd, ed, l):
	
	# Defining variables
	c_info = c.strip().split("_")
	circ_name = ("%s_5linseq" % c)
	ENST = c_info[3]
	exon2_rank = ex[0]			# Rank of circRNA exon
	exon1_rank = ex[0] - 1		# Rank of exon before circRNA
	exon2_pos = exon2_rank -1	# Position in sorted_transcript_dic of circRNA exon
	exon1_pos = exon1_rank -1	# Position in sorted_transcript_dic of exon before circRNA
	exon2_ENSE = rd[ENST][exon2_pos][1]
	exon1_ENSE = rd[ENST][exon1_pos][1]
	
	exon2_searchterm = str("%s_%s" % (ENST, exon2_ENSE))	# Search term for circRNA exon in exon_dic
	exon1_searchterm = str("%s_%s" % (ENST, exon1_ENSE))	# Search term for exon before the circRNA in exon_dic
		
	# Searchin exon information
	exon2_info = ed[exon2_searchterm]
	exon1_info = ed[exon1_searchterm]
	
	
	# Defining coordinates of exons
	exon2_start = exon2_info[4]
	exon2_end = exon2_info[5]
	exon1_start = exon1_info[4]
	exon1_end = exon1_info[5]
	strand = ex[2]
	
	final_output = [circ_name, l, exon1_ENSE, exon1_start, exon1_end, exon2_ENSE, exon2_start, exon2_end, strand]
	print("%s" % ('\t'.join(final_output)))

"""

"""
FUNCTION: seq3
OBJECTIVE: obtain the coordinates of the last exon of the circRNA and
the coordinates of the exon after it.

c: refers to circ
ex: refers to exon (obtained from sorted_circ_dic[circ][0])
rd: refers to sorted_rank_dic
ed: refers to exon_dic
l: refers to label
"""

"""
def seq3(c, ex, rd, ed, l):
	
	# Defining variables
	c_info = c.strip().split("_")
	circ_name = ("%s_3linseq" % c)
	ENST = c_info[3]
	exon1_rank = ex[0]			# Rank of circRNA exon
	exon2_rank = ex[0] + 1		# Rank of after before circRNA
	exon1_pos = exon1_rank -1	# Position in sorted_transcript_dic of circRNA exon
	exon2_pos = exon2_rank -1	# Position in sorted_transcript_dic of exon after circRNA
	exon1_ENSE = rd[ENST][exon1_pos][1]
	exon2_ENSE = rd[ENST][exon2_pos][1]
	
	exon1_searchterm = str("%s_%s" % (ENST, exon1_ENSE))	# Search term for circRNA exon in exon_dic
	exon2_searchterm = str("%s_%s" % (ENST, exon2_ENSE))	# Search term for exon after the circRNA in exon_dic
		
	# Searchin exon information
	exon1_info = ed[exon1_searchterm]
	exon2_info = ed[exon2_searchterm]
	
	
	# Defining coordinates of exons
	exon1_start = exon1_info[4]
	exon1_end = exon1_info[5]
	exon2_start = exon2_info[4]
	exon2_end = exon2_info[5]
	strand = ex[2]
	
	final_output = [circ_name, l, exon1_ENSE, exon1_start, exon1_end, exon2_ENSE, exon2_start, exon2_end, strand]
	print("%s" % ('\t'.join(final_output)))
"""

########################################################################
########################################################################
########################################################################
########################################################################

# Main cycle for exons
for exon in exon_data:
	data = exon.strip().split("\t")
	save_exons(data, exon_dic) # Fills exon_dic
	save_ranks(data, rank_dic) # Fills rank_dic


# Main cycle for circRNA overlaps
for circ in overlap_data:
	data = circ.strip().split("\t")
	circ_transcripts(data, circ_dic, unaligned_exon_dic, intertranscript_dic, exon_dic)

# Sorting circ_dic
for circ in circ_dic:
	ex_list = circ_dic[circ]
	
	def getKey(item):
		return item[0]
	
	# Filling sorted dictionary
	sorted_circ_dic[circ] = sorted(ex_list, key=getKey)

# Purging unsorted circRNA dictionary
circ_dic = {}

# Sorting rank_dic
for transcript in rank_dic:
	rank_list = rank_dic[transcript]
	
	def getKey(item):
		return item[0]
	
	# Filling sorted dictionary
	sorted_rank_dic[transcript] = sorted(rank_list, key=getKey)

# Purging unsorted rank dictionary
rank_dic = {}

"""
# Checks for linear sequences to report
for circ in sorted_circ_dic:
	
	# Defining transcript
	transcript = circ.strip().split("_")[3]
	
	# Defining min and max rank of transcript
	minimum = min(sorted_rank_dic[transcript])
	maximum = max(sorted_rank_dic[transcript])
	rank_min = minimum[0]
	rank_max = maximum[0]
	
	####################################################################
	# Check for multiexonic circRNAs
	if len(sorted_circ_dic[circ]) == 2:
		
		# Calculating ranks of circ exons
		rank_ex5 = sorted_circ_dic[circ][0][0]
		rank_ex3 = sorted_circ_dic[circ][1][0]
		
		# Checking position of circRNA in the transcript
		if rank_ex5 > rank_min and rank_ex3 < rank_max:

			label = "multiexon"
			seq5(circ, sorted_circ_dic[circ][0], sorted_rank_dic, exon_dic, label)
			seq3(circ, sorted_circ_dic[circ][1], sorted_rank_dic, exon_dic, label)
		
		elif rank_ex5 == rank_min and rank_ex3 < rank_max:

			label = "multiexon_5boundary" # Multiexon at 5' of transcript
			seq3(circ, sorted_circ_dic[circ][1], sorted_rank_dic, exon_dic, label)
		
		elif rank_ex5 > rank_min and rank_ex3 == rank_max:

			label = "multiexon_3boundary" # Multiexon at 3' of transcript
			seq5(circ, sorted_circ_dic[circ][0], sorted_rank_dic, exon_dic, label)
			
		elif rank_ex5 == rank_min and rank_ex3 == rank_max:
			label = "missing_exon"
			print("%s\t%s" % (circ, label))

	####################################################################
	# Check for monoexonic circRNAs
	elif len(sorted_circ_dic[circ]) == 1:
		
		# Calculating rank of circ exon
		rank = sorted_circ_dic[circ][0][0]

		# Checks for true monoexonic circRNA
		circ_start = int(circ.strip().split("_")[1])
		circ_end = int(circ.strip().split("_")[2])
		ENSE = sorted_circ_dic[circ][0][1]
		ENST = circ.strip().split("_")[3]
		exon_search = ("%s_%s" % (ENST, ENSE))
		
		exon_start = int(exon_dic[exon_search][4])
		exon_end = int(exon_dic[exon_search][5])
		
		if circ_start != exon_start or circ_end != exon_end:
			label = "missing_exon"
			print("%s\t%s" % (circ, label))
		
		else:	
		
			if rank > rank_min and rank < rank_max:
				# Monoexonic circ with both sequences calculated
				label = "monoexon"
				seq5(circ, sorted_circ_dic[circ][0], sorted_rank_dic, exon_dic, label)
				seq3(circ, sorted_circ_dic[circ][0], sorted_rank_dic, exon_dic, label)
			
			elif rank == rank_min and rank < rank_max:
				# Monoexon at 5' of transcript
				label = "monoexon_5boundary" 
				seq3(circ, sorted_circ_dic[circ][0], sorted_rank_dic, exon_dic, label)
			
			elif rank > rank_min and rank == rank_max:
				# Monoexon at 3' of transcript
				label = "monoexon_3boundary"
				seq5(circ, sorted_circ_dic[circ][0], sorted_rank_dic, exon_dic, label)
	
	####################################################################
	# Check for missing exons
	else:
		label = "missing_exon"
		print("%s\t%s" % (circ, label))
"""

####################### CIRCULARIZING EXONS SECTION ####################


for circ in sorted_circ_dic:
	if len(sorted_circ_dic[circ]) == 2:
		
		rank1 = sorted_circ_dic[circ][0][0]
		rank2 = sorted_circ_dic[circ][1][0]
		
		if (len(rank1) == 3 and len(rank2) == 3):
			label = "multiexon"
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
			(circ, label,
			
			sorted_circ_dic[circ][0][3], # ENSG of first exon
			sorted_circ_dic[circ][0][4], # ENST of first exon	
			sorted_circ_dic[circ][0][1], # ENSE of first exon
			sorted_circ_dic[circ][0][0], # Rank of first exon
			sorted_circ_dic[circ][0][2], # Strand of first exon
			
			sorted_circ_dic[circ][1][3], # ENSG of second exon
			sorted_circ_dic[circ][1][4], # ENST of second exon	
			sorted_circ_dic[circ][1][1], # ENSE of second exon
			sorted_circ_dic[circ][1][0], # Rank of second exon
			sorted_circ_dic[circ][1][2])) # Strand of second exon
		
		elif (len(rank1) == 4 or len(rank2) == 4):
			
			label = "intronic"
			print("%s\t%s" % (circ, label))
	
	elif len(sorted_circ_dic[circ]) == 1:
		
		# Calculating rank of circ exon
		rank = sorted_circ_dic[circ][0][0]

		# Checks for true monoexonic circRNA
		circ_start = int(circ.strip().split("_")[1])
		circ_end = int(circ.strip().split("_")[2])
		ENSE = sorted_circ_dic[circ][0][1]
		ENST = circ.strip().split("_")[3]
		exon_search = ("%s_%s" % (ENST, ENSE))
		
		exon_start = int(exon_dic[exon_search][4])
		exon_end = int(exon_dic[exon_search][5])
		
		if circ_start != exon_start or circ_end != exon_end:
			
			# Fetching gene start and gene end
			gene_start = int(exon_dic[exon_search][8])
			gene_end = int(exon_dic[exon_search][9])
			
			# Checking for intronic or intergenic circRNAs
			if (circ_start >= gene_start and circ_end <= gene_end):
				
				if (circ_start < exon_start and len(rank) == 3):
					
					if circ in unaligned_exon_dic:
						
						label=("putativeexon")
						
						# Check needed to write the output in the correct order
						if (unaligned_exon_dic[circ][0][2] == '+'):
						
							print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
							(circ, label,
							
							unaligned_exon_dic[circ][0][3], # ENSG of first exon
							unaligned_exon_dic[circ][0][4], # ENST of first exon	
							unaligned_exon_dic[circ][0][1], # ENSE of first exon
							unaligned_exon_dic[circ][0][0], # Rank of first exon
							unaligned_exon_dic[circ][0][2], # Strand of first exon
							
							sorted_circ_dic[circ][0][3], # ENSG of second exon
							sorted_circ_dic[circ][0][4], # ENST of second exon	
							sorted_circ_dic[circ][0][1], # ENSE of second exon
							sorted_circ_dic[circ][0][0], # Rank of second exon
							sorted_circ_dic[circ][0][2])) # Strand of second exon
						
						else:
							print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
							(circ, label,
							
							sorted_circ_dic[circ][0][3], # ENSG of second exon
							sorted_circ_dic[circ][0][4], # ENST of second exon	
							sorted_circ_dic[circ][0][1], # ENSE of second exon
							sorted_circ_dic[circ][0][0], # Rank of second exon
							sorted_circ_dic[circ][0][2], # Strand of second exon
							
							unaligned_exon_dic[circ][0][3], # ENSG of first exon
							unaligned_exon_dic[circ][0][4], # ENST of first exon	
							unaligned_exon_dic[circ][0][1], # ENSE of first exon
							unaligned_exon_dic[circ][0][0], # Rank of first exon
							unaligned_exon_dic[circ][0][2])) # Strand of first exon
							
					else:
						label = "intronic"
						print("%s\t%s" % (circ, label))
					
				
				elif (circ_end > exon_end and len(rank) == 3):
					
					if circ in unaligned_exon_dic:
						
						label=("putativeexon")
						
						# Check needed to write the output in the correct order
						if (unaligned_exon_dic[circ][0][2] == '+'):
							
							print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
							(circ, label,
	
							sorted_circ_dic[circ][0][3], # ENSG of second exon
							sorted_circ_dic[circ][0][4], # ENST of second exon	
							sorted_circ_dic[circ][0][1], # ENSE of second exon
							sorted_circ_dic[circ][0][0], # Rank of second exon
							sorted_circ_dic[circ][0][2], # Strand of second exon
							
							unaligned_exon_dic[circ][0][3], # ENSG of first exon
							unaligned_exon_dic[circ][0][4], # ENST of first exon	
							unaligned_exon_dic[circ][0][1], # ENSE of first exon
							unaligned_exon_dic[circ][0][0], # Rank of first exon
							unaligned_exon_dic[circ][0][2])) # Strand of first exon
						
						else:
							
							print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
							(circ, label,
								
							unaligned_exon_dic[circ][0][3], # ENSG of first exon
							unaligned_exon_dic[circ][0][4], # ENST of first exon	
							unaligned_exon_dic[circ][0][1], # ENSE of first exon
							unaligned_exon_dic[circ][0][0], # Rank of first exon
							unaligned_exon_dic[circ][0][2], # Strand of first exon
							
							sorted_circ_dic[circ][0][3], # ENSG of second exon
							sorted_circ_dic[circ][0][4], # ENST of second exon	
							sorted_circ_dic[circ][0][1], # ENSE of second exon
							sorted_circ_dic[circ][0][0], # Rank of second exon
							sorted_circ_dic[circ][0][2])) # Strand of second exon
					
					else:
						label = "intronic"
						print("%s\t%s" % (circ, label))
				
				
				else:
					label = "intronic"
					print("%s\t%s" % (circ, label))
							
			
			elif (circ_start < gene_start or circ_end > gene_end):
				label = "intergenic"
				print("%s\t%s" % (circ, label))
		
		elif len(rank) == 3:
			label = "monoexon"
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (circ,
			label,
			
			sorted_circ_dic[circ][0][3], # ENSG of first exon
			sorted_circ_dic[circ][0][4], # ENST of first exon
			sorted_circ_dic[circ][0][1], # ENSE of exon
			sorted_circ_dic[circ][0][0], # Rank of exon
			sorted_circ_dic[circ][0][2])) # Strand of exon
		
		elif len(rank) == 4:
			label = "intronic"
			print("%s\t%s" % (circ, label))
	
	elif len(sorted_circ_dic[circ]) == 0:
		
		# Defining variables needed to check for intronic or intergenic condition
		circ_start = int(circ.strip().split("_")[1])
		circ_end = int(circ.strip().split("_")[2])

		ENST = circ.strip().split("_")[3]
	
		# A complete key needs to be reconstructed in order to search
		# the exon dictionary correctly. Any exon of the ENST considered
		# will be ok, since the searched gene values are all the same
		# for every exon in the same transcript
		
		ENSE = sorted_rank_dic[ENST][0][1]
		exon_search = "%s_%s" % (ENST, ENSE)
		
		"""
		# Old method for INSANE FUCKS THAT DO NOT REMEMBER THE FUCKING
		# DICTIONARIES THAT THEY MADE THEMSELVES
		
		for key in exon_dic.keys():
			if key.startswith(ENST):
				exon_search = key
				break
		"""

		gene_start = int(exon_dic[exon_search][8])
		gene_end = int(exon_dic[exon_search][9])
		transcript_start = int(exon_dic[exon_search][10])
		transcript_end = int(exon_dic[exon_search][11])

		# Checking for intronic or intergenic circRNA
		if (circ_start >= gene_start and circ_end <= gene_end):
			
			# NEWCODE
			
			if (circ in unaligned_exon_dic and circ_start >= transcript_start and circ_end <= transcript_end):
				label = "putativeexon"
				print("%s\t%s\t%s\t%s\t%s\t%s\t%s" %
				(circ, label,

				unaligned_exon_dic[circ][0][3], # ENSG of first exon
				unaligned_exon_dic[circ][0][4], # ENST of first exon	
				unaligned_exon_dic[circ][0][1], # ENSE of first exon
				unaligned_exon_dic[circ][0][0], # Rank of first exon
				unaligned_exon_dic[circ][0][2])) # Strand of first exon
				
			
			else:
				label = "intronic"
				print("%s\t%s" % (circ, label))
		
		elif (circ_start < gene_start or circ_end > gene_end):
						
			label = "intergenic"
			print("%s\t%s" % (circ, label))
	
	else:
		print("ERROR - len of sorted_circ_dic is %s" % (len(sorted_circ_dic[circ])))
	
# Printing intertranscript circRNAs
for circ in intertranscript_dic:
	if len(intertranscript_dic[circ]) == 2:
		print("%s_%s\tintergenic" % (circ, intertranscript_dic[circ][0][1]))
		print("%s_%s\tintergenic" % (circ, intertranscript_dic[circ][1][1]))
                
                # It was decided to refer to intertranscript circs as intergenic
                # The old printing statements are below
		# print("%s_%s\tintertranscript" % (circ, intertranscript_dic[circ][0][1]))
		# print("%s_%s\tintertranscript" % (circ, intertranscript_dic[circ][1][1]))
        
################ DEBUG SECTION ########################################

"""
print("\nexon_dic")
for exon in exon_dic:
	print("%s\t%s" % (exon, exon_dic[exon]))

	
print("\nsorted_circ_dic")
for circ in sorted_circ_dic:
	print("%s\t%s" % (circ, sorted_circ_dic[circ]))

print("\nsorted_rank_dic")
for rank in sorted_rank_dic:
	print("%s\t%s" % (rank, sorted_rank_dic[rank]))

print("\nunaligned_exon_dic")
for exon in unaligned_exon_dic:
	print("%s\t%s" % (exon, unaligned_exon_dic[exon]))

print("\nintertranscript_dic")
for circ in intertranscript_dic:
	print("%s\t%s" % (circ, intertranscript_dic[circ]))
"""
