#!/usr/bin/python3

from Bio import SeqIO
import argparse
import enum
import os, sys
import subprocess
#import csv

#nimblegen create-project -p PROJECT_FOLDER -s SAMPLE_FOLDER -g REFERENCE_FOLDER --regexp_paired_end
#nimblegen align -p PROJECT_FOLDER <params>
#nimblegen metrics -p PROJECT_FOLDER
#nimblegen methyl-analysis -p PROJECT_FOLDER -c CONTROL_GENOME


#docker run -v `pwd`/case-study:/data -u `id -u \`whoami\`` -it pytest:latest create-project -p progetto -s data/samples -r data/ref/lambda.fa --regexp _R1.fastq

#alias nimblegen="docker run -v `pwd`/case-study:/data -u `id -u \`whoami\`` -it pynimble"
#nimblegen create-project -p progetto -s data/samples -r data/ref/lambda.fa --regexp _R1.fastq
#nimblegen  methyl-analysis -p progetto -c J02459.1

#TODO - se progetto esiste gia, lancia FileExistsError -- gestire!!

class script(enum.Enum):
    alignment = "bsmap.sh"
    filtering = "filtering.sh"
    rmduplicates = "rmdups.sh"
    metrics = "metrics.sh"
    methratio = "ch3_analysis.sh"
    bissnp = "bissnp.sh"

class folders(enum.Enum):
    alignment = "alignments/"
    workspace = "workspace/"
    metrics = "metrics/"
    methylation = "methylation/"


class inputs(enum.Enum):
    sample = "sample"
    reference = "ref"
    primary_target = "primary"
    capture_target = "capture"
    dbsnp = "dbsnp"


class Nimblegen(object):
    def __init__(self, project):
        self.project = project
        self.reference = None
        self.primary = None
        self.capture = None
        self.dbsnp = None
        self.samples = dict()

        self.project_folder = "/data/{}".format(self.project)
        self.files_project = "{}/project.txt".format(self.project_folder)
        self.command_log = "{}/commands.txt".format(self.project_folder) #TODO -> da usare
        self.output_path = "{}/output".format(self.project_folder)
        self.path_shcript = "/bin/" #bash script  (?)

    def chromosome_sizes(self):
        chromosomes = dict()

        for fasta_record in SeqIO.parse(self.reference, "fasta"):
            chromosomes[fasta_record.id] = len(fasta_record)

        return chromosomes


    def create_project(self, samples_folder, regexp, reference, primary_target, capture_target, dbsnp):
        os.makedirs(self.project_folder)

        samples = self.preprocess_samples(samples_folder, regexp)
#        chr_sizes = self.chromosome_sizes() #TODO - salvare da qualche altra parte

        # crea file di configurazione
        self.write_project(samples, inputs.sample)
        self.write_project(reference, inputs.reference)
        self.write_project(primary_target, inputs.primary_target)
        self.write_project(capture_target, inputs.capture_target)
        self.write_project(dbsnp, inputs.dbsnp)

        # crea directory output e relative sottodirectory
        os.makedirs(self.output_path)

        for out in folders:
            os.makedirs("{}/{}".format(self.output_path, out.value))


    def alignment(self):
        self.read_project()

        print("INFO: performing read alignment against the reference genome")

        alignment_folder = "{}/{}".format(self.output_path, folders.alignment.value)
        workspace_folder = "{}/{}".format(self.output_path, folders.workspace.value)
        metrics_folder = "{}/{}".format(self.output_path, folders.metrics.value)

        for sample, data in self.samples.items():
            r1 = "{}/{}".format(data["path"], data["r1"])
            r2 = "{}/{}".format(data["path"], data["r2"])

            params = [sample, r1, r2, self.reference, alignment_folder]

            if self.exec_script(script.alignment, params) != 0:
                print("Something is happened during alignment of sample {}".format(sample))
        else:
            print("allineamento ok")

        #rimozione duplicati
        print("INFO: performing duplicate removal")
        params = [alignment_folder, workspace_folder, metrics_folder]

        if self.exec_script(script.rmduplicates, params) != 0:
            print("CROSTICINA RMDUPS")

        print("INFO: performing filtering for mapped and properly paired-end reads + clipping overlapping mates")
        params = [workspace_folder, alignment_folder]

        if self.exec_script(script.filtering, params) != 0:
            print("KRAVARINCI")

    def metrics_assessment(self):
        self.read_project()

        alignment_folder = "{}/{}".format(self.output_path, folders.alignment.value)
        workspace_folder = "{}/{}".format(self.output_path, folders.workspace.value)
        metrics_folder = "{}/{}".format(self.output_path, folders.metrics.value)

        print("INFO: assessing metrics")
        params = [alignment_folder, self.reference, workspace_folder, metrics_folder, self.primary, self.capture]

        if self.exec_script(script.metrics, params) != 0:
            print("E TI PAREVA. PD")


    def methylation_extraction(self, control_genome):
        self.read_project()

        alignment_folder = "{}/{}".format(self.output_path, folders.alignment.value)
        workspace_folder = "{}/{}".format(self.output_path, folders.workspace.value)
        metrics_folder = "{}/{}".format(self.output_path, folders.metrics.value)
        methylation_folder = "{}/{}".format(self.output_path, folders.methylation.value)

        print("\nINFO: methylation level and bisulfite conversion efficiency estimations\n")
        params = [alignment_folder, self.reference, workspace_folder, control_genome, methylation_folder]

        if self.exec_script(script.methratio, params) != 0:
            print("ESPLOSO TUTTO")


        print("\nINFO: methylation calling using BisSNP tool\n")
        params = [alignment_folder, self.reference, workspace_folder, methylation_folder, self.dbsnp, self.capture]

        if self.exec_script(script.bissnp, params) != 0:
            print("CRISTO IN CROCE")


    def exec_script(self, which, params):
        path_script = "/bin/./{}".format(which.value)

        return subprocess.call([path_script] + params)

    def read_project(self):
        self.samples = dict()

        with open(self.files_project) as f:
            for line in f:
                fields = line.strip().split("\t")
                tag = fields[0]

                if tag == inputs.sample.value:
                    sample, r1, r2, path = fields[1:]

                    self.samples[sample] = {
                        "path": path,
                        "r1": r1,
                        "r2": r2
                    }
                elif tag == inputs.reference.value:
                    self.reference = fields[1]
                elif tag == inputs.primary_target.value:
                    self.primary = fields[1]
                elif tag == inputs.capture_target.value:
                    self.capture = fields[1]
                elif tag == inputs.dbsnp.value:
                    self.dbsnp = fields[1]


    def write_project(self, what, tag):
        tag = tag.value

        with open(self.files_project, "a") as f:
            if type(what) is str:
                f.write("{}\t/{}\n".format(tag, what))

            if type(what) is dict:
                for sample, data in what.items():
                    r1, r2 = data["r1"], data["r2"]
                    path = data["path"]

                    f.write("{}\t{}\t{}\t{}\t/{}\n".format(tag, sample, r1, r2, path))

    def preprocess_samples(self, folder, regexp):
        """ Returns samples information about R1 and R2 fastq file """
        samples = dict()

        for sample_dir in os.listdir(folder):  # for each directory
            current_path = "{}/{}".format(folder, sample_dir)
            sample = os.listdir("{}".format(current_path))

            r1_sample = next(filter(lambda x: regexp in x, sample))
            samplename = r1_sample.split(regexp)[0]

            r2_sample = next(filter(lambda x: regexp not in x, sample))

            samples[samplename] = {
                "path": current_path,
                "r1": r1_sample,
                "r2": r2_sample
            }

        return samples



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="nimblegen-pipeline")
    subparsers = parser.add_subparsers(help="...")

    parser_project = subparsers.add_parser("create-project", help="create-project help")
    parser_alignment = subparsers.add_parser("align", help="alignment help")
    parser_metrics = subparsers.add_parser("metrics", help="alignment help")
    parser_analysis = subparsers.add_parser("methyl-analysis", help="alignment help")

    parser_project.add_argument("-p", "--project", dest="project_folder", metavar="project-folder",
                        action="store", type=str, required=True)
    parser_project.add_argument("-s", "--samples", dest="samples_folder", metavar="samples-folder",
                        action="store", type=str, required=True)
    parser_project.add_argument("-r", "--reference", dest="reference_folder", metavar="reference-folder",
                        action="store", type=str, required=True)
    parser_project.add_argument("--dbsnp", dest="dbsnp", metavar="dbsnp-vcf-file",
                        action="store", type=str, required=True)
    parser_project.add_argument("--primary", dest="primary_target", metavar="primary-target",
                        action="store", type=str, required=True)
    parser_project.add_argument("--capture", dest="capture_target", metavar="capture-target",
                        action="store", type=str, required=True)
    parser_project.add_argument("--regexp", dest="regexp_paired_end", action="store", type=str)


    parser_alignment.add_argument("-p", "--project", dest="project_folder", metavar="project-folder",
                        action="store", type=str, required=True)

    parser_metrics.add_argument("-p", "--project", dest="project_folder", metavar="project-folder",
                        action="store", type=str, required=True)

    parser_analysis.add_argument("-p", "--project", dest="project_folder", metavar="project-folder",
                        action="store", type=str, required=True)

    parser_analysis.add_argument("-c", "--control_genome", dest="control_genome", metavar="control-genome",
                        action="store", type=str, required=True)

    args = parser.parse_args()

    try:
        command = sys.argv[1]
    except IndexError:
        parser.print_help(sys.stderr)
        sys.exit()

    nimble = Nimblegen(args.project_folder)

    if command == "create-project":
        nimble.create_project(args.samples_folder, \
            args.regexp_paired_end, args.reference_folder, \
            args.primary_target, args.capture_target, \
            args.dbsnp)

    elif command == "align":
        nimble.alignment()

    elif command == "metrics":
        nimble.metrics_assessment()

    elif command == "methyl-analysis":
        nimble.methylation_extraction(args.control_genome)
