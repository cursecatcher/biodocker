#!/usr/bin/env python

import os,sys

"""
DATA EXTRACTOR
Extracts the data of file 1 requested from file 2
"""

# Defining data dictionary
data_dic ={}

# Opening input files
data_file = open(sys.argv[1]).readlines()
query_file = open(sys.argv[2]).readlines()

# Defining index
data_index = int(sys.argv[3])
query_index = int(sys.argv[4])


"""
Function: fill_dic
Fills data_dictionary using a user-defined column as source for keys
"""
def fill_dic (x, index, dic):

	dkey = x[index] # Defining dictionary key using user supplied value
	
	if dkey not in dic:
		dic[dkey] = x

"""
Function: search_data
Obtains data from data_dic based on a user-defined query index
"""
def search_data (x, q, dic):
	
	query_key = x[q] # Querykey
	
	if query_key in dic:
		out = dic[query_key] + x # Building output

	else:
		repeats = int(len(dic[dic.keys()[0]]))
		out = ["NA"]*repeats + x
	
	print("\t".join(out))

# Main cycle
for line in data_file:

	data = line.strip().split("\t")
	fill_dic(data, data_index, data_dic) # Filling data dictionary

for query_line in query_file:

	query_data = query_line.strip().split("\t")
	search_data(query_data, query_index, data_dic) # Obtaining data

# DEBUG
"""
print("Printing data file")
for line in data_file:
	data = line.strip().split("\t")
	print(data)

print(data_index)


for info in data_dic:
	print("%s\t%s" % (info, data_dic[info]))
"""
