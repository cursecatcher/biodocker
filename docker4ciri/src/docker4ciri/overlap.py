#!/usr/bin/env python3

# This script calculates the overlap between two or more files containing circRNA data

import argparse
import csv 
import os 


def get_tool_from_filename(filename):
    """ Extract the name of the tool from the filename of the file produced by reformat script"""
    tokens = os.path.basename(filename).split(".")
    return (filename, tokens[-2].lower()) if len(tokens) >= 3 \
            else (filename, None)

def get_boolean_row(expected_header, used_tools):
    """ Return a True/False list. The i-th element represent the presence 
    of the i-th tool of expected_header in used_tools """
    return [str(tool in used_tools)[0] for tool in expected_header]


def fix_coordinates(coordinates, used_tool):
    if used_tool in ("circexplorer2",):
        chromo, start, end, circ_id, strand = coordinates
        start = str(int(start) + 1)
        circ_id = "{}_{}_{}".format(chromo, start, end)
        coordinates = chromo, start, end, circ_id, strand 
        
    return coordinates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", dest="input_files", action="store", nargs="*", type=str, required=True)
    parser.add_argument("-o", "--out", dest="output_dir", action="store", type=str, required=True)
    parser.add_argument("-t", "--threshold", dest="min_support", action="store", default=1, type=int)

    args = parser.parse_args() 

    input_files = [get_tool_from_filename(filename) for filename in args.input_files]

    #collecting tools from input filenames
    tools = set(tool for filename, tool in input_files if tool is not None)

    circdict = dict() 
    #collect circRNAs 
    for input_file, used_tool in input_files: 
        if used_tool is None:
            continue 

        print("Processing {}".format(input_file))

        with open(input_file) as f: 
            for line in csv.reader(f, delimiter="\t"): 
                mycirc = fix_coordinates(tuple(line), used_tool)
                
                if mycirc not in circdict:
                    circdict[mycirc] = {used_tool}
                else:
                    circdict[mycirc].add(used_tool)

    #filter circRNAs detected by a sufficient number of tools 
    with open(os.path.join(args.output_dir, "circRNA_detection.txt"), "w") as overlap_file,\
         open(os.path.join(args.output_dir, "circRNA_list.txt"), "w") as circrna_list_file:

        sorted_tools = list(tools)
        
        #for each circRNA, how many tools (and which) detect it (?)
        overlap_csv = csv.writer(overlap_file, delimiter="\t")
        overlap_csv.writerow(["ID", "#detection"] + sorted_tools)
        #circRNA list detected at least by N tools 
        csv_circrna_list = csv.writer(circrna_list_file, delimiter="\t") #no header 
        rows = list() 

        for circ, tools in circdict.items():
            #keeping circRNAs detected at least by n tools (summary)
            if len(tools) >= args.min_support:
                csv_circrna_list.writerow(circ)

            boolean_flags = get_boolean_row(sorted_tools, tools)  
            circ_id = "{}_{}".format(circ[3], "+" if circ[-1] == "1" else "-")
            n_evidences = boolean_flags.count("T")
            
            rows.append([circ_id, n_evidences] + boolean_flags)
        
        #sorting circ data by #detection 
        rows.sort(key=lambda l: l[1], reverse=True)
        overlap_csv.writerows(rows)
        