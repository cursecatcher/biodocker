#!/usr/bin/env python3
import argparse
from datetime import datetime 
import enum 
import os 
import subprocess
import sys 

class KnownPaths(enum.Enum):
    inputFolder = "/data"
    alignmentFolder = "/data/{}_alignments"
    referenceFolder = "/ref"
    scratchFolder = "/scratch"

class Sample(object):
    def __init__(self, samplename, r1, r2):
        self.name = samplename 
        self.r1_fastq = r1 
        self.r2_fastq = r2 

def get_trimmed_name(fastq_file):
    i = fastq_file.index(".")
    return fastq_file[:i] + ".trimmed" + fastq_file[i:]

def get_samplename(fastq_file, regexp):
    index = fastq_file.index(regexp)
    return fastq_file[:index]    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    parser_bsmap = subparser.add_parser("bsmap")
    parser_bismark = subparser.add_parser("bismark")
    parser_bsmap.add_argument("-rgxR1", dest="regexp_R1", type=str, required=True)
    parser_bismark.add_argument("-rgxR1", dest="regexp_R1", type=str, required=True)

    args, unknownargs = parser.parse_known_args()
    unknownargs_str = " ".join(unknownargs)

    try:
        mode = sys.argv[1]
        regexp_R1 = args.regexp_R1
    except IndexError:
        parser.print_help()
        #show bsmap / bismark -h 
        sys.exit(2)
    
    fastq_files = set(
        fastq for fastq in os.listdir(KnownPaths.inputFolder.value) \
            if fastq.endswith(".fastq") or fastq.endswith(".fastq.gz")
    )
    if len(fastq_files) == 0:
        sys.exit("No FASTQ file found in mounted folder!")
    r1_files = set(fastq_name for fastq_name in fastq_files if regexp_R1 in fastq_name)
    r2_files = fastq_files - r1_files

    samples = list() 
    
    for fastq in r1_files:
        prefix_R1 =  get_samplename(fastq, regexp_R1)
        mate = None 
        try:
            mate = [e for e in r2_files if prefix_R1 in e][0]
            print("Found mate of {} -> {}\nProducing output file {}".format(fastq, mate, get_trimmed_name(mate)), flush=True)
        except IndexError:
            print("Cannot find a mate for {} file".format(fastq), flush=True)
        if mate is not None: 
            fastq_r1 = os.path.join(KnownPaths.inputFolder.value, fastq)
            fastq_r2 = os.path.join(KnownPaths.inputFolder.value, mate)
            sample = Sample(prefix_R1, fastq_r1, fastq_r2)
            samples.append(sample)
    
    alignmentFolder = KnownPaths.alignmentFolder.value.format(mode)
    try:
        os.makedirs(alignmentFolder)
    except FileExistsError:
        print("Folder {} already exists. Overriding data...".format(alignmentFolder), flush=True)

    if mode == "bsmap":
        try:
            fasta = [f for f in os.listdir(KnownPaths.referenceFolder.value) if f.endswith(".fa") or f.endswith(".fna")]
            if len(fasta) > 1:
                sys.exit("More than one FASTA file present in the mounted folder. Unsupported.")
            fasta = fasta.pop()
        except FileNotFoundError:
            sys.exit("Cannot locate reference file. Please ensure you've mount FASTA directory in /ref")
        except IndexError:
            sys.exit("No FASTA file found in mounted directory.")

        fasta = os.path.join(KnownPaths.referenceFolder.value, fasta)
        command = "bsmap -a {} -b {} -d {} -o {} {}"

        for sample in samples: 
            samfile = os.path.join(alignmentFolder, sample.name + ".sam")
            current_command = command.format(
                sample.r1_fastq, sample.r2_fastq, fasta, samfile, unknownargs_str
            )
            subprocess.run(current_command, shell=True)
            
    elif mode == "bismark":
        if not os.path.exists(KnownPaths.scratchFolder.value):
            sys.exit("Bismark requires a temporary folder to store temporary files.\nPlease mount a scratch folder in /scratch before run this docker.")

        command = "bismark --bowtie2 --genome_folder {} --temp {} -1 {} -2 {} --output_dir {} {}"
        ref = KnownPaths.referenceFolder.value
        flag = True 

        while flag:
            try:
                temp_folder = os.path.join(
                    KnownPaths.scratchFolder.value, 
                    mode + "_" + datetime.now().strftime("%d_%m_%Y__%H_%M_%S")
                )
                os.makedirs(temp_folder)
                flag = False 
            except FileExistsError:
                pass 

        for sample in samples:
            current_command = command.format(
                ref, temp_folder, sample.r1_fastq, sample.r2_fastq, alignmentFolder, unknownargs_str
            )
            subprocess.run(current_command, shell=True)
        
    subprocess.run("chmod -R 777 {}".format(alignmentFolder), shell=True)