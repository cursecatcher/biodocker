#!/usr/bin/env python3

import argparse
import csv 
import enum 
import os, sys 
from circutils import * 


class ToolParameters(object):
    def __init__(self):
        #TODO - complete headers 
        self.__parameters = {
            SupportedTool.ACFS:
                (None, 
                (0, 1, 2, 5)), 
            SupportedTool.CIRI:
                (["circRNA_ID", "chr", "circRNA_start", "circRNA_end", "#junction_reads", 
                "SM_MS_SMS", "#non_junction_reads", "junction_reads_ratio", "circRNA_type", 
                "gene_id", "junction_reads_ID"], 
                (1, 2, 3, -1)), 
            SupportedTool.CIRI2:
                (["circRNA_ID", "chr", "circRNA_start", "circRNA_end", "#junction_reads", "SM_MS_SMS",
                "#non_junction_reads", "junction_reads_ratio",	"circRNA_type", "gene_id", "strand", "junction_reads_ID"], 
                (1, 2, 3, 10)), 
            SupportedTool.CIRCEXPLORER:
                (["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", 
                "itemRgb", "exonCount", "exonSizes", "exonOffsets", "readNumber", "circType", 
                "geneName", "isoformName", "exonIndex/intronIndex", "flankIntron"], 
                (0, 1, 2, 5)), 
            SupportedTool.CIRCEXPLORER2:
                (None, 
                (0, 1, 2, 5)), 
            SupportedTool.STARCHIP: 
                (["cRNA","P1GeneInfo","P1Strand","P1GeneDistance","P2GeneInfo", "P2Strand", "P2GeneDistance"], 
                (0, -1, -1, 2)), #cRNA field contains chromosome, start and end information
            SupportedTool.UROBORUS:
                (None,  
                (0, 1, 2, 3)), 
            SupportedTool.CIRCRNAFINDER: 
                (["chr", "start", "end", "name", "score", "strand"], 
				(0, 1, 2, 5)),
            SupportedTool.FINDCIRC2:
                (["chrom", "start", "end", "name", "n_reads", "strand", "n_uniq", "uniq_bridges", "best_qual_left",
				"best_qual_right ", "tissue", "tiss_counts", "edits", "anchor_overlap", "breakpoints", "signal", "strandmatch", "category"],  
				(0, 1, 2, 5)),
            SupportedTool.KNIFE:
                #field #0 contains the whole circRNA coordinate: chr15|RPS17:82823288|RPS17L:83207122|rev|- 
                (["junction", "orig_count", "orig_posterior", "swapped_count", "swapped_posterior", "total_reads"],  
				(0,-1, -1, -1)),
            SupportedTool.DCC: 
                (["chr", "start", "end", "genename", "junctiontype", "strand", "circRNA region", "overall regions"],
				(0,1,2,5))   
        }

    def __getitem__(self, key):
        return self.__parameters[key]


# This script reformats the output files produced by various circRNA detector tools
# in order to obtain a common format for further analysis. 
# The output file will be saved in the same directory of the input file. 
# The output file extension will contain the information regarding the tool used to detect the circRNAs. 

def check_chromosome(chr_value):
    """Reformat the chromosome value. Return a string in the format chrN  """

    chr_value = str(chr_value)#.lower()

    if "chr" not in chr_value:
        chr_value = f"chr{chr_value}"
    
    return chr_value

def check_strand(strand_value):
    """Reformat the strand value. Return a string representing 1 if the strand is the forward one, else -1. """
    strand = str(strand_value) #forcing string type 

    if strand in ("+", "-"):
        strand = "1" if strand == "+" else "-1"

    return strand

def process_line(line, indexes, tool):
    """ Return a 4-tuple containing the info about chromosome, start and end indexes and strand info, if it present """

    if tool == SupportedTool.KNIFE:
        #special case #1: all the info are in the same field 
        #e.g. : chr13|DNAJC3:96377506|DNAJC3:96375496|rev|+
        return process_knife_line(line[indexes[0]])

    elif tool == SupportedTool.STARCHIP:
        #special case #2: chromosome info and start/end indexes are in the same field,
        #while the strand info has a dedicated field
        #e.g. : chr10:103622470-103661114
        chromo, coordinates = line[indexes[0]].split(":")
        start, end = coordinates.split("-")
        strand = line[indexes[-1]]
        return chromo, start, end, strand 
    else:
        return tuple(line[index] if index >= 0 else 1 for index in indexes)

def reformat(filename, expected_header, indexes, tool):
    """ This function is a generator: it reformats the input file one line at the time, 
    yielding the formatted version of current line. """ 

    with open(filename) as fi:
        circsv = csv.reader(fi, delimiter="\t")   
        actual_header = next(circsv, None)

        if expected_header is not None and actual_header != expected_header:
            sys.exit("The actual header does not match the expected one. Check if you have selected the right tool and/or input file.")
        else: 
            #there is no header - return to the beginning of file 
            fi.seek(0, 0)

        for line in circsv:
            #extract fields from line 
            chr, begin, end, strand = process_line(line, indexes, tool)
            #standardize some stuff 
            chr = check_chromosome(chr) 
            strand = check_strand(strand) 
            #save result 
            yield [chr, begin, end, "{}_{}_{}".format(chr, begin, end), strand]


# input file: output file produced by one of the supported tools 
# output file: filename of the input file reformatted in a common format

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", dest="input_path", action="store", required=True)
    parser.add_argument("-o", "--out", dest="output_path", action="store", required=True)
    parser.add_argument("-t", "--tool", dest="used_tool", action="store", required=True, choices=[tool.value for tool in SupportedTool])

    args = parser.parse_args()
    used_tool = SupportedTool.get_enum(args.used_tool)

    if used_tool not in SupportedTool:
        sys.exit("Seems that the tool {} is not supported yet.".format(args.used_tool))

    for current_file in os.listdir(args.input_path):
        print("Reformatting {}".format(current_file))

        current_file_path = os.path.join(args.input_path, current_file)
        filenamepath, extension = os.path.splitext(current_file_path)
        filename = os.path.basename(filenamepath)
        output_file = os.path.join(args.output_path, "{}.{}.reformat".format(filename, args.used_tool))

        with open(output_file, "w") as fo: 
            csvw = csv.writer(fo, delimiter="\t")
            header, indexes = ToolParameters()[used_tool]

            for line in reformat(current_file_path, header, indexes, used_tool):
                csvw.writerow(line)


        

