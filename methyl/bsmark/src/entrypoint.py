#!/usr/bin/env python3 

import enum 
import subprocess
import sys 

class Tools(enum.Enum):
    bismark = "bismark"
    bsmap = "bsmap"

def find_bam_and_sam_files(args):
    return [arg for arg in args if arg.endswith(".bam") or arg.endswith(".sam")]

def get_options(args):
    return [arg for arg in args if arg.startswith("-")]


if __name__ == "__main__":
    my_tool = sys.argv[1]

    if my_tool in [tool.value for tool in Tools]:
        command_line_args = " ".join(sys.argv[1:])
        command = "python3 /usr/bin/bsmark.py {}".format(command_line_args)
        print("Running aligment with command: {}".format(command))
        subprocess.call(command, shell=True)
        sys.exit(0)
    

    if my_tool == "deduplicate_bismark":
        #fai samtools sort , poi deduplicate_bismark 
        sorted_bams = list() 
        options = get_options(sys.argv[1:])

        if "-h" in options:
            subprocess.run("deduplicate_bismark -h", shell=True, check=True)
            sys.exit(0)

        for a_file in find_bam_and_sam_files(sys.argv[1:]):
            bam_file = a_file.replace(".sam", ".bam") if a_file.endswith(".sam") else a_file
            sorted_file = "sorted_{}".format(bam_file)

            print("Sorting {} in temporary file {}...".format(a_file, sorted_file))
            command = "samtools sort -n {} > {}".format(a_file, sorted_file)
            ret_code = subprocess.run(command, shell=True, check=True)

            if ret_code.returncode != 0:
                print("Something gone wrong during sorting of {}".format(a_file))
                continue

            sorted_bams.append(sorted_file)
        
        if len(sorted_bams) > 0:
            bam_files = " ".join(sorted_bams)
            options = " ".join(options)
            command = "deduplicate_bismark {} {}".format(options, bam_files)

            print("Launching the following command: {}".format(command))
            subprocess.run(command, check=True, shell=True)
        else:
            print("No sorted BAM found in the current directory")

    
    else:
        command = " ".join(sys.argv[1:])
        print("Launching command {}".format(command))
        subprocess.run(command, check=True, shell=True)
    