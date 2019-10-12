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
                (None, (0, 1, 2, 5)), 
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
                (None, (0, 1, 2, 5)), 
            # SupportedTool.STARCHIP: ([], ()), 
            SupportedTool.UROBORUS:
                (None,  (0, 1, 2, 3))
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
        chr_value = "chr" + chr_value
    
    return chr_value

def check_strand(strand_value):
    """Reformat the strand value. Return a string representing 1 if the strand is the forward one, else -1. """
    strand = str(strand_value) #forcing string type 

    if strand in ("+", "-"):
        strand = "1" if strand == "+" else "-1"

    return strand


def reformat(filename, expected_header, indexes):
    """ This function is a generator: it reformats the input file one line at the time, 
    yielding the formatted version of current line. """ 

    chr_index, start_index, end_index, strand_index = indexes 

    with open(filename) as fi:
        circsv = csv.reader(fi, delimiter="\t")   
        actual_header = next(circsv, None)

        if expected_header is not None and actual_header != expected_header:
            sys.exit("The actual header does not match the expected one. Check if you have selected the right tool.")

        indexes = [chr_index, start_index, end_index, strand_index] 

        for line in circsv:
            #extract fields from line 
            chr, begin, end, strand = [line[index] if index >= 0 else 1 for index in indexes]
            #standardize some stuff 
            chr = check_chromosome(chr) 
            strand = check_strand(strand) 
            #save result 
            yield [chr, begin, end, "{}_{}_{}".format(chr, begin, end), strand]


# input file: output file produced by one of the supported tools 
# output file: filename of the input file reformatted in a common format

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", dest="input_file", action="store", required=True)
    parser.add_argument("-o", "--out", dest="output_path", action="store", required=True)
    parser.add_argument("-t", "--tool", dest="used_tool", action="store", required=True, choices=[tool.value for tool in SupportedTool])

    args = parser.parse_args()
    used_tool = SupportedTool.get_enum(args.used_tool)

    if used_tool not in SupportedTool:
        sys.exit("Seems that the tool {} is not supported yet.".format(args.used_tool))

    filenamepath, extension = os.path.splitext(args.input_file)
    filename = os.path.basename(filenamepath)
    output_file = os.path.join(args.output_path, "{}.{}.circrna".format(filename, args.used_tool))

    with open(output_file, "w") as fo: 
        csvw = csv.writer(fo, delimiter="\t")
        header, indexes = ToolParameters()[used_tool]

        for line in reformat(args.input_file, header, indexes):
            csvw.writerow(line)
    

