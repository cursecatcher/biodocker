#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os, os.path as path
import circRNAannot
import enum


class AnnotationSources(enum.Enum):
    CIRCBASE = ("circbase", circRNAannot.CircBaseDB)
    TSCD = ("tscd", circRNAannot.tscdDB)
    EXORBASE = ("exorbase", circRNAannot.ExoRBaseDB)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in", "-i", dest="circrna", action="store", type=str, required=True)
    parser.add_argument("--scratch", "-s", dest="scratch_folder", action="store", type=str, required=True)
    parser.add_argument("--out", "-o", dest="output_folder", action="store", type=str, default=None)
    parser.add_argument("--ref", "-r", dest="assembly_version", action="store", type=str, choices=[x.value for x in circRNAannot.AssemblyVersion], required=True)
    parser.add_argument("--sources", dest="sources", action="store", nargs="*", type=str, choices=[x.value[0] for x in AnnotationSources])

    args = parser.parse_args()
    
    if args.sources is None: 
        args.sources = [x.value[0] for x in AnnotationSources]

    output_folder = path.dirname(path.abspath(args.circrna)) if args.output_folder is None \
                    else path.abspath(args.output_folder)

    circset = circRNAannot.CircRNAs()
#    circset.add_circ(circRNA("", "chr6","4891946", "4892613", "+"))
    assembly = args.assembly_version

    with open(args.circrna) as f:
        print("Loading circRNAs from {}".format(args.circrna))
        csvr = csv.reader(f, delimiter="\t")

        for line in csvr:
            chromo, start, end, _, strand = line
            strand = "+" if int(strand) == 1 else "-"
            circset.add_circ(chromo, start, end, strand)

    for source in AnnotationSources: 
        db_name, db_class = source.value 

        print("Annotating circRNAs from {}".format(db_name))
        db_class(args.assembly_version)\
            .download_db(args.scratch_folder)\
                .annotate(circset)\
                    .save(output_folder)

    # if AnnotationSources.CIRCBASE.value in args.sources:
    #     circbase = circRNAannot.CircBaseDB(args.assembly_version)
    #     circbase.download_db(args.scratch_folder)
    #     circbase.annotate(circset).save(output_folder)

    # if AnnotationSources.TSCD.value in args.sources:
    #     for version in ("adult", "fetal"):
    #         tscd = circRNAannot.tscdDB(args.assembly_version, version)
    #         tscd.download_db(args.scratch_folder)
    #         tscd.annotate(circset).save(output_folder)
        
    # if AnnotationSources.EXORBASE.value in args.sources: 
    #     exorbase = circRNAannot.ExoRBaseDB(args.assembly_version)
    #     exorbase.download_db(args.scratch_folder)
    #     exorbase.annotate(circset).save(output_folder)


    #delete scratch folder content
    if False:
        scratch = path.abspath(args.scratch_folder)
        for filename in os.listdir(scratch):
            trash = path.join(scratch, filename)

            if path.exists(trash):
                print("Deleting {}".format(filename))
                os.remove(trash)
