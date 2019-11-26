#!/usr/bin/env python3

import argparse
import csv 
import os 

# This script calculates the overlap between two or more files containing circRNA data

def get_tool_from_filename(filename):
    """ Extract the name of the tool from the filename of the file produced by reformat script"""
    return os.path.basename(filename).split(".")[-2] 

def get_boolean_row(expected_header, used_tools):
    """ Return a True/False list. The i-th element represent the presence 
    of the i-th tool of expected_header in used_tools """
    return [str(tool in used_tools)[0] for tool in expected_header]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", dest="input_files", action="store", nargs="*", type=str, required=True)
    parser.add_argument("-o", "--out", dest="output_dir", action="store", type=str, required=True)
    parser.add_argument("-t", "--threshold", dest="min_support", action="store", default=1, type=int)

    args = parser.parse_args() 

    #collecting tools from input filenames
    tools = set(get_tool_from_filename(filename) for filename in args.input_files)

    circdict = dict() 
    #collect circRNAs 
    for input_file in args.input_files: 
        print("Loading {}".format(input_file))
        used_tool = get_tool_from_filename(input_file)

        with open(input_file) as f: 
            for line in csv.reader(f, delimiter="\t"): 
                mycirc = tuple(line) 
                
                if mycirc not in circdict:
                    circdict[mycirc] = set() 

                circdict[mycirc].add(used_tool)

    #filter circRNAs detected by a sufficient number of tools 
    with open(os.path.join(args.output_dir, "circRNA_detection.txt"), "w") as summary,\
         open(os.path.join(args.output_dir, "circRNA_list.txt"), "w") as detected:

        sorted_tools = list(tools)
        print(sorted_tools)
        
        #for each circRNA, how many tools (and which) detect it (?)
        summary_csv = csv.writer(summary, delimiter="\t")
        summary_csv.writerow(["ID", "#detection"] + sorted_tools)
        #circRNA list detected at least by N tools 
        detected_csv = csv.writer(detected, delimiter="\t") #no header 
        rows = list() 

        for circ, tools in circdict.items():
            #keeping circRNAs detected at least by n tools (summary)
            if len(tools) >= args.min_support:
                detected_csv.writerow(circ)

            boolean_flags = get_boolean_row(sorted_tools, tools)            
            circ_id = "{}_{}".format(circ[3], "+" if circ[-1] == 1 else "-")
            n_evidences = boolean_flags.count("T")
            
            rows.append([circ_id, n_evidences] + boolean_flags)
        
        #sorting circ data by #detection 
        rows.sort(key=lambda l: l[1], reverse=True)
        summary_csv.writerows(rows)
        