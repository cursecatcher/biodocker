#!/usr/bin/env python3

import argparse
import datetime
import enum
import os
import re
import shutil
import subprocess
import sys

class Commands(enum.Enum):
    merge_data = "merge_data"
    merge = "merge"
    annotation = "annotation"
    structure = "structure"
    ciri2 = "ciri2"


class Scripts(enum.Enum):
    ciri_as = "/ciri2/CIRI_AS_v1.2.pl"
    ciri2 = "/ciri2/CIRI2.pl"
    annotation = "/ciri2/annotate_circrna.py"
    ciri_merge = "/ciri2/merge_v2.py"
    data_merge = "/ciri2/deseq_merge.py"

class Paths(enum.Enum):
    #main scratch folder: a subfolder will be created for each execution
    scratch_folder = "/scratch"
    #path where output files will be stored
    output_folder = "/data/out"
    #input file for CIRI_AS and annotation procedures
    cirifile = "/data/cirifile"
    #input file for CIRI_AS and CIRI2 procedures
    samfile = "/data/samfile"
    #input file in FASTA format for CIRI_AS and CIRI2 procedures
    reference = "/data/reference"
    #annotation file for CIRI_AS and CIRI2 procedures
    annotation = "/data/annotation"
    #input folder containing samples to merge (ciri_merge and merge_data procedures)
    folder_to_merge = "/data/input_merge"

class OutputFiles(enum.Enum):
    ciri_merge = "merged_circRNA"
    data_merge = "MergedData.tsv"
    ciri2_prediction = "prediction.ciri"
    ciri_as_prediction = "structure"
#    annotation = "chiaro"

def getAnnotationExtension(path):
    tokens = path.split("/")
    data_folder = "/{}".format("/".join(tokens[:-1]))

    filename = [f for f in os.listdir(data_folder) if tokens[-1] in f].pop()
    return filename.split(".")[-1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    merge_data = subparsers.add_parser(Commands.merge_data.value)
    merge_data.add_argument("-s", "--samples", dest="samples", nargs="*", required=True)
    merge_data.add_argument("--cov", nargs="*", dest="covariates", required=True)
    merge_data.add_argument("--order", dest="cov_order", nargs="*")
    merge_data.add_argument("--col", dest="column", type=int, default=1)
    merge_data.add_argument("--ext", dest="extension", default="ciri")

    merge = subparsers.add_parser(Commands.merge.value)
    merge.add_argument("-s", "--samples", dest="samples", nargs="*", required=True)
    merge.add_argument("--cov", nargs="*", dest="covariates", required=True)
    merge.add_argument("--order", dest="cov_order", nargs="*")
    merge.add_argument("--mr", dest="min_reads", type=int, required=True)
    merge.add_argument("--mrep", dest="min_replicates", type=int, required=True)
    merge.add_argument("--avg", dest="min_average", type=int, required=True)

    annotation = subparsers.add_parser(Commands.annotation.value)
    annotation.add_argument("-s", "--sources", dest="annotation_sources", nargs="*", required=True)
    annotation.add_argument("-v", "--version", dest="assembly_version", type=str, default="hg19")

    ciri_as = subparsers.add_parser(Commands.structure.value)
    ciri_as.add_argument("-a", "--anno", dest="annotation", action="store_true")

    ciri2 = subparsers.add_parser(Commands.ciri2.value)
    ciri2.add_argument("-S", "--max_span", dest="max_span", type=int, default=200000)
    ciri2.add_argument("-T", "--thread_num", dest="nthreads", type=int, default=1)
    ciri2.add_argument("-U", "--mapq_uni", dest="quality_threshold", type=int, default=10)
    ciri2.add_argument("--strigency", dest="strigency", choices=["high", "low", "zero"], default="high")
    ciri2.add_argument("-a", "--anno", dest="annotation", action="store_true")

    args = parser.parse_args()

    mode = sys.argv[1] if len(sys.argv) > 1 else None

    if mode is None:
        parser.print_help(sys.stderr)
        print("No mode specified. Bye.")
        sys.exit(1)

    #################################################
    #local scratch folder creation
    #################################################
    scratch_id = re.sub("[:\.\ ]", "-", str(datetime.datetime.now()))
    my_scratch = "{}/{}".format(Paths.scratch_folder.value, scratch_id)

    try:
        os.makedirs(my_scratch)
    except OSError:
        print("Error: cannot create scratch directory {}. Try again.".format(scratch_id))
        sys.exit()
        #TODO: fix

    print("Current scratch folder ID: {}".format(scratch_id), flush=True)


    if mode == Commands.merge.value:
        #set default order
        if args.cov_order is None:
            args.cov_order = sorted(list(set(args.covariates)))
        #do merge
        command = [
            "python3", Scripts.ciri_merge.value,
            "--in {}".format(Paths.output_folder.value),
#            "--in {}".format(Paths.folder_to_merge.value),
            "--out {}/{}".format(Paths.output_folder.value, OutputFiles.ciri_merge.value),
            "--samples {}".format(" ".join(args.samples)),
            "--cov {}".format(" ".join(args.covariates)),
            "--order {}".format(" ".join(args.cov_order)),
            "--mr {}".format(args.min_reads),
            "--mrep {}".format(args.min_replicates),
            "--avg {}".format(args.min_average)
        ]

    elif mode == Commands.merge_data.value:
        #set default order
        if args.cov_order is None:
            args.cov_order = sorted(list(set(args.covariates)))
        #do merge
        command = [
            "python3", Scripts.data_merge.value,
            "--in {}".format(Paths.output_folder.value),
#            "--in {}".format(Paths.folder_to_merge.value),
            "--out {}/{}".format(Paths.output_folder.value, OutputFiles.data_merge.value),
            "--samples {}".format(" ".join(args.samples)),
            "--cov {}".format(" ".join(args.covariates)),
            "--order {}".format(" ".join(args.cov_order)),
            "--col {}".format(args.column),
            "--ext {}".format(args.extension)
        ]

    elif mode == Commands.annotation.value:
        command = [
            "python3", Scripts.annotation.value,
            "--in {}".format(Paths.cirifile.value),
            "--out {}".format(Paths.output_folder.value),
            "--scratch {}".format(my_scratch),
            "--sources {}".format(" ".join(args.annotation_sources)),
            "--ref {}".format(args.assembly_version)
        ]

    elif mode == Commands.structure.value:
        command = [
            "perl", Scripts.ciri_as.value,
            "--sam {}".format(Paths.samfile.value),
            "--ciri {}".format(Paths.cirifile.value),
            "--ref_file {}".format(Paths.reference.value),
            "--out {}/{}".format(Paths.output_folder.value, OutputFiles.ciri_as_prediction.value)
        ]
        if args.annotation:
            extension = getAnnotationExtension(Paths.annotation.value)
            command.append("--anno {}.{}".format(Paths.annotation.value, extension))

    elif mode == Commands.ciri2.value:
        strigency_value = "--{}_strigency".format(
            "no" if args.strigency == "zero" else args.strigency
        )

        command = [
            "perl", Scripts.ciri2.value,
            "--in {}".format(Paths.samfile.value),
            "--out {}/{}".format(Paths.output_folder.value, OutputFiles.ciri2_prediction.value),
            "--ref_file {}".format(Paths.reference.value),
            "--max_span {}".format(args.max_span),
            "--thread_num {}".format(args.nthreads),
            "--mapq_uni {}".format(args.quality_threshold),
            strigency_value
        ]
        if args.annotation:
            extension = getAnnotationExtension(Paths.annotation.value)
            command.append("--anno {}.{}".format(Paths.annotation.value, extension))

    #execute task
    ret = subprocess.run(" ".join(command), shell=True)

    #remove scratch folder
    if ret.returncode == 0:
        
        print("Removing scratch folder")
        shutil.rmtree(my_scratch)

    sys.exit(ret.returncode)