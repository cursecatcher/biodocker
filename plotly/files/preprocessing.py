#!/usr/bin/env python3

import csv
import sys
from math import log10

#heatmap.sh nomefile
## crea scratch
## chiama preprocessing.py : scratch, separator, status
### crea count table filtrata e formattata in scratch
## chiama heatmap.R: scratch,
### produce heatmap.html in scratch, sposta in outputfolder
## cancella scratch

def get_genename(genename):
    return genename.strip().replace("\"", "")
#    return genename.split(":")[1] if ":" in genename else genename

def logify(num_list):
    l = None
    try:
        l = [log10(float(x)) for x in num_list]
    except ValueError:
        print("Error: trying to obtain log10 of a value <= 0!! Check your input file!", file=sys.stderr)
    return l


if __name__ == "__main__":
    inputfile = "/data/in/count_table"
    genelist = "/data/in/gene_list"
    scratchfolder = sys.argv[1]
    separator = sys.argv[2].upper()
    raw_counts = sys.argv[3].upper() == "RAW"

    filteredfile = "{}/myfile.csv".format(scratchfolder)

    if separator == "TAB":
        my_separator = "\t"
    elif separator == "COMMA":
        my_separator = ","
    elif separator == "SPACE":
        my_separator = " "
    else:
        print("Unsupported field separator.", file=sys.stderr)
        sys.exit(1)

    geneset = set()

    print("Filtering count table...")

    with open(genelist) as genefile:
        for line in genefile:
            geneset.add(get_genename(line))

    print("Producing filtered file...")

    with open(inputfile) as fi, open(filteredfile, "w") as fo:
        csvr = csv.reader(fi, delimiter=my_separator)
        csvw = csv.writer(fo, delimiter="\t")

        #write header
        csvw.writerow(next(csvr, None))

        for line in csvr:
            genename = get_genename(line[0])
            genedata = line[1:]

            if genename in geneset:
                if raw_counts:
                    genedata = logify(genedata)

                    if genedata is None:
                        print("Cannot create filtered count table because of incompatible values. Abort", file=sys.stderr)
                        sys.exit(2)

                csvw.writerow([genename] + list(genedata))

    print("File {} created".format(filteredfile))

    sys.exit(0)
