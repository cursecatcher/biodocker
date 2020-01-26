#!/usr/bin/env python2.7

##UNIVOCAL CLASSIFIER
##This script outputs a univocal circRNA classification

import os
import sys

def dic_debugger (dic):
	"""
	Function: dic_debugger
	"""

	for entry in dic:
		print("%s\t%s" % (entry, dic[entry]))


def class_converter (x):
	"""
	Function: class_converter
	Converts circRNA classification in a sortable string

	Every classification beginning with 01 will have priority over
	the main isoform if the main isoform itself does not contain it.
	"""

	converter = {
		"multiexon":        "01_multiexon",
		"monoexon":         "01_monoexon",
		"putativeexon":     "01_putativeexon",
		"intronic":         "05_intronic",
		"intergenic":       "06_intergenic",
		"intertranscript":  "07_intertranscript",
	}

	return converter.get(x)



def data_parser (x):
	"""
	Function: data_parser
	Organize data and pass it to dic_filler function
	"""

	y = x[2].split("_") # Splitting circRNA name
	z = x[0].split("-") # Splitting gene name
	iso_index = x[1]
	ENST = y[3]

	# Variables for dic_circ (univocal circRNA dic)
	#circ_name = "%s_%s_%s" % (y[0], y[1], y[2]) # Defining circRNA name
	circ_name = "%s_%s_%s_%s" % (y[0], y[1], y[2], z[0]) # Defining circRNA name
	class_index = class_converter(x[3])
	circ_values = (iso_index, class_index, ENST)

	# Variables for circRNA raw data inside dic_data
	line_name = "%s_%s" % (x[2], x[3])
	line_data = x

	return (circ_name, circ_values, line_name, line_data)


def fill_dic(dic, k, v, check_presence):
	"""
	Fills a dictionary with provided key and values
	check_presence	= True: check if the key is already in the dictionary
					= False: doesn't check dictionary (use only when data append is not needed)
	"""

	if check_presence:
		# Key existence check
		# Data is appended to the key
		if k not in dic:
			dic[k] = []
		dic[k].append(v)
	else:
		# No check for key existence
		# Data overwrite if duplicate keys are present
		dic[k] = v


def dic_sorter (dic, dic_sorted):
	"""
	Function: dic_sorter
	Sorts the tuples of the various keys in a dictionary
	"""

	for key in dic:
		dic_sorted[key] = sorted(dic[key])


def classifier (dic, key, udic):
	"""
	Function: circRNA_classifier
	Allows to obtain a univocal classification for a circRNA
	"""

	if key not in udic:
		udic[key] = []

	# First pass for higher importance classifications
	for element in dic[key]:
		class_index = element[1].strip().split("_")

		# Appending to univocal dictionary only if it is empty
		if str(class_index[0]) == "01" and len(udic[key]) == 0:
			udic[key] = element

	# Second pass for lower importance classifications
	for element in dic[key]:

		# Appending to univocal dictionary onlu if no higher importance classification was found
		if len(udic[key]) == 0 and str(class_index[0]) != "07":
			udic[key] = element

	# Second pass for intertranscript classifications
	for element in dic[key]:

		# Appending to univocal dictionary onlu if no higher importance classification was found
		if len(udic[key]) == 0:
			udic[key] = element


def circRNA_data_retrival (dic, key):
	"""
	Function: circRNA_data_retrival
	Returns detailed information about the selected circRNA isoform
	"""

	# Reconstructing key for info dictionary
	cn = key.strip().split("_")							# Splitting for circ name
	circ_name = ("%s_%s_%s" % (cn[0], cn[1], cn[2]))	# Obtaining chr_start_end
	ENST = dic[key][2]									# Obtaining ENST
	classification = dic[key][1].strip().split("_")[1]	# Obtaining classification

	circ_key = ("%s_%s_%s" % (circ_name, ENST, classification))

	return (circ_key)


def output_parser (dic, key):
	"""
	Function: output parser
	"""

	circRNA_name = ["_".join(circ_key.strip().split("_")[0:3])]	# circRNA name (chr_start_end)
	classification = [circ_data[3]]								# Classification
	isoform_info = circ_data[0:2]								# Transcript name
	other_data = circ_data[4:]

	output_data = circRNA_name + classification + isoform_info + other_data
	return (output_data)


# MAIN PROGRAM
if __name__ == "__main__":
	classification_dic = {}			#circRNA classification dictionary
	sorted_classification_dic = {}	#circRNA classification dictionary SORTED
	data_dic = {}					#circRNA data retrival dictionary
	univocal_dic = {}				#Final classification dictionary

	# Opening input file containing CircHunter's circRNA classification
	with open(sys.argv[1]) as circ_file:
		for circ in circ_file:
			data = circ.strip().split("\t")
			dparsed = data_parser(data) # Extract data from line
			fill_dic(classification_dic, dparsed[0], dparsed[1], True) #Fill classification dictionary
			fill_dic(data_dic, dparsed[2], dparsed[3], False) # Fill data dictionary with circRNA full line data

	dic_sorter(classification_dic, sorted_classification_dic) # Sorting classification dic
	#classification_dic = {} # Cleaning up memory


	# Univocal classification
	for circRNA in sorted_classification_dic:
		classifier(sorted_classification_dic, circRNA, univocal_dic)
	#sorted_classification_dic = {} # Cleaning up memory

	# Data retrival
	for circRNA in univocal_dic:
		circ_key = circRNA_data_retrival(univocal_dic, circRNA)	# Retriving circRNA key
		circ_data = data_dic[circ_key]							# Retriving circRNA data
		final_data = output_parser (circ_data, circ_key)		# Formatting output
		print("\t".join(final_data))


	# Debug
	# print("\nDEBUG")
	# print("\nClassification dic")
	# dic_debugger(classification_dic)
	# print("\nSORTED Classification dic")
	# dic_debugger(sorted_classification_dic)
	# print("\nData dic")
	# dic_debugger(data_dic)
	# print("\nUnivocal dic")
	# dic_debugger(univocal_dic)
