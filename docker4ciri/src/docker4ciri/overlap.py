#!/usr/bin/env python3

import argparse
import csv 
import os 

# This script calculates the overlap between two or more files containing circRNA stuff
#TODO - continuare questa bellissima descrizione 

def get_tool_from_filename(filename):
    """ Extract the name of the tool from the filename of the file produced by reformat script"""
    return os.path.basename(filename).split(".")[-2] 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", dest="input_files", action="store", nargs="*", type=str, required=True)
    parser.add_argument("-t", "--threshold", dest="min_support", action="store", default=0, type=int)

    args = parser.parse_args() 

    circdict = dict() 

    #collect circRNAs 
    for input_file in args.input_files: 
        used_tool = get_tool_from_filename(input_file)

        with open(input_file) as f: 
            for line in csv.reader(f, delimiter="\t"): 
                
                mycirc = tuple(line) 
                if mycirc not in circdict:
                    circdict[mycirc] = list() 
                circdict[mycirc].append(used_tool)

    #filter circRNAs detected by a sufficient number of tools 
    for circ, tools in circdict.items():
        if len(tools) >= args.min_support:
            tools = [",".join(tools)]
            #print(" ".join(list(circ) + tools))
            #TODO - scrivi su file magari -> in che formato? 
            pass 
 
