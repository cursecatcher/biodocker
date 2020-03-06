#!/usr/bin/env python3
import argparse
import enum 
import os 
import subprocess
import sys 


class KnownPaths(enum.Enum):
    inputFolder = "/data"
    trimFolder = "/data/trimmed"
    
    #trimFolder = "/data/trimmed"

def get_trimmed_name(fastq_file):
    i = fastq_file.index(".")
    return fastq_file[:i] + ".trimmed" + fastq_file[i:]

def get_samplename(fastq_file, regexp):
    index = fastq_file.index(regexp)
    return fastq_file[:index]    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    parser_qc = subparser.add_parser("qc")
    parser_cut = subparser.add_parser("cutadapt")
    parser_cut.add_argument("-rgxR1", dest="regexp_R1", type=str, required=True)
    #parse input/output options in order to ignore them ? 

    args, unknownargs = parser.parse_known_args()

    try:
        mode = sys.argv[1]
    except IndexError:
        parser.print_help()
        #show cutadapt -h / qc -h 
        sys.exit(2)
    
    fastq_files = set(
        fastq for fastq in os.listdir(KnownPaths.inputFolder.value) \
            if fastq.endswith(".fastq") or fastq.endswith(".fastq.gz")
    )

    if len(fastq_files) == 0:
        sys.exit("No FASTQ file found in mounted folder!")

    if mode == "qc":
        #fastqc +  multiqc 
        command = "start_qc.sh {}".format(" ".join(sys.argv[2:]))
        subprocess.run(command, shell=True)

    elif mode == "cutadapt":
        regexp = args.regexp_R1
        r1_fastq = set(fastq_name for fastq_name in fastq_files if regexp in fastq_name)
        r2_fastq = fastq_files - r1_fastq

        #pairing mates 
        for fastq in r1_fastq:
            prefix_R1 =  get_samplename(fastq, regexp)
            mate = None 
            try:
                mate = [e for e in r2_fastq if prefix_R1 in e][0]
                print("Found mate of {} -> {}".format(fastq, mate))
                print(get_trimmed_name(mate))
            except IndexError:
                print("Cannot find a mate for {} file".format(fastq), flush=True)
            if mate is not None: 
                fastq_r1 = os.path.join(KnownPaths.inputFolder.value, fastq)
                fastq_r2 = os.path.join(KnownPaths.inputFolder.value, mate)
                trimmed_r1 = os.path.join(KnownPaths.trimFolder.value, get_trimmed_name(fastq))
                trimmed_r2 = os.path.join(KnownPaths.trimFolder.value, get_trimmed_name(mate))

                cutadapt_args = " ".join(unknownargs)
                
                os.makedirs(KnownPaths.trimFolder.value, exist_ok=True)
                subprocess.run("chmod 777 {}".format(KnownPaths.trimFolder.value), shell=True)

                command = f"cutadapt {cutadapt_args} -o {trimmed_r1} -p {trimmed_r2} {fastq_r1} {fastq_r2}"
                print("Running the following command: {}".format(command), flush=True)
                subprocess.run(command, shell=True)
            

    