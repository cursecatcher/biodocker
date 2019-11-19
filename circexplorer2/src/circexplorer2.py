#!/usr/bin/env python2 

import argparse
import os, sys 
import enum 
import subprocess
import shlex

class Commands(enum.Enum):
    Align = "align"
    Annotate = "annotate"
    Parse = "parse"
    # Assemble = "assemble"
    # Denovo = "denovo"

class Folders(enum.Enum):
    input_dir = "/data/in"
    output_dir = "/data/out"

class MountPoints(object):
    genome_fasta = os.path.join(Folders.input_dir.value, "fasta") 
    index_bowtie1 = os.path.join(Folders.input_dir.value, "index_bowtie1") 
    index_bowtie2 = os.path.join(Folders.input_dir.value, "index_bowtie2") 
    gtf_file = os.path.join(Folders.input_dir.value, "gtf_file")
    fastq_file = os.path.join(Folders.input_dir.value, "fastq")
    junctions_file = os.path.join(Folders.input_dir.value, "junctions")
    annotation_file = os.path.join(Folders.input_dir.value, "annotation")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    aligner_c = subparsers.add_parser(Commands.Align.value)
    annotator_c = subparsers.add_parser(Commands.Annotate.value)
    parser_c = subparsers.add_parser(Commands.Parse.value)
    
    aligner_c.add_argument("-p", dest="thread_number", action="store", type=int, default=10)
    aligner_c.add_argument("-o", dest="output_filename", action="store", type=str)
    aligner_c.add_argument("-r", dest="reference_type", action="store", type=str, choices=["fasta", "bowtie1", "bowtie2"], required=True)

    annotator_c.add_argument("-o", dest="output_filename", action="store", type=str)

    parser_c.add_argument("-t", dest="aligner", action="store", type=str, choices=["BWA", "STAR"], required=True)
    parser_c.add_argument("-b", "--bed", dest="output_filename", action="store", type=str, default="back_spliced_junction.bed")

    args = parser.parse_args()

    try: 
        mode = sys.argv[1]
    except IndexError:
        parser.print_help()
        sys.exit(1)

    if mode == Commands.Align.value:
        if args.reference_type == "fasta":
            reference = "-g {}".format(MountPoints.genome_fasta)
        elif args.reference_type == "bowtie1":
            reference = "-i {}".format(MountPoints.index_bowtie1)
        elif args.reference_type == "bowtie2":
            reference = "-j {}".format(MountPoints.index_bowtie2)
        
        command_tail = [
            "{}".format(reference),
            "-G {}".format(MountPoints.gtf_file), 
            "-f {}".format(MountPoints.fastq_file), 
            "-p {}".format(args.thread_number)
        ]

    elif mode == Commands.Annotate.value: 
        command_tail = [
            "-r {}".format(MountPoints.annotation_file), 
            "-g {}".format(MountPoints.genome_fasta), 
            "-b {}".format(MountPoints.junctions_file), 
            "-o {}".format(args.output_filename)
        ]
    
    elif mode == Commands.Parse.value: 
        command_tail = [
            "-t {}".format(args.aligner), 
            "-b {}".format(os.path.join(Folders.output_dir.value, args.output_filename)),
            MountPoints.junctions_file
        ]

    command_tail = " ".join(command_tail)
    logfile = os.path.join(Folders.output_dir.value, "log_{}.txt".format(mode))
    command = "CIRCexplorer2 {} {} > {}".format(mode, command_tail, logfile)

    print("Running CIRCexplorer2 docker in {} mode. Running the following command:\n{}\n".format(mode, command))
    
    ret = subprocess.call(command, shell=True)
    print("Return value: {}".format(ret))

    if ret == 0:
        print("CIRCexplorer2 terminated successfully")
    else: 
        subprocess.call("CIRCexplorer2 {} -h".format(mode), shell=True)
        sys.exit(1)

    