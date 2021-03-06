#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os, os.path as path
import circRNAannot
import enum
import traceback
import urllib

class AnnotationSources(enum.Enum):
    CIRCBASE = ("circbase", circRNAannot.CircBaseDB)
    TSCD = ("tscd", circRNAannot.tscdDB)
    EXORBASE = ("exorbase", circRNAannot.ExoRBaseDB)
    CIRCRIC = ("circric", circRNAannot.CircRicDB)
    CSCD = ("cscd", circRNAannot.CSCDDB)
    CIRC2DISEASE = ("circ2disease", circRNAannot.Circ2DiseaseDB)
    CIRCFUNBASE = ("circfunbase", circRNAannot.CircFunBaseDB)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in", "-i", dest="circrna", action="store", type=str, required=True)
    parser.add_argument("--scratch", "-s", dest="scratch_folder", action="store", type=str, required=True)
    parser.add_argument("--out", "-o", dest="output_folder", action="store", type=str, default=None)
    parser.add_argument("--ref", "-r", dest="assembly_version", action="store", type=str, 
        choices=[x.value for x in circRNAannot.AssemblyVersion], required=True)
    parser.add_argument("--sources", dest="sources", action="store", nargs="*", type=str, 
        choices=[x.value[0] for x in AnnotationSources])

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
        
        if db_name in args.sources:
            print("Annotating data from {} db".format(db_name))
            try:
                db_class(assembly)\
                    .download_db(args.scratch_folder)\
                        .annotate(circset)\
                            .save(output_folder)

            except circRNAannot.NoDatabaseException: 
                assembly = circRNAannot.AssemblyVersion.get_enum_value(assembly)
                organism = circRNAannot.Organism.get_organism_by_assembly(assembly)
                print("\033[01;31mCannot annotate from {} db: it is not available for {} organism.\033[00m".format(db_name, organism))
            except urllib.error.HTTPError: 
                print("Cannot download {} db because it is unreachable. Please try later.".format(db_name))
            except Exception as e:
                print("Some very random event happened: {}".format(e)) #TODO - provvisorio, fix 
                print("###########################")
                traceback.print_exc()
                print("###########################")


    #delete scratch folder content
    if False:
        scratch = path.abspath(args.scratch_folder)
        for filename in os.listdir(scratch):
            trash = path.join(scratch, filename)

            if path.exists(trash):
                print("Deleting {}".format(filename))
                os.remove(trash)
