#!/usr/bin/env python3 

import argparse
import csv
import os, sys
from circutils import process_knife_line

# This script merges 2+ CIRI2/CIRCexplorer2 output files. 
# It produces the following output files: 
# Count table: rows = circRNAs, columns = samples 
# CircRNA file: report all the predicted circRNAs 


class sample_group(object):
    """Class sample_group describes a group of samples in the same condition"""

    def __init__(self, group_id):
        self.__id = group_id
        self.__samples = list()

    @property
    def samples(self):
        return list(self.__samples)

    @property
    def group_id(self):
        return self.__id

    def add_sample(self, sample_id):
        self.__samples.append(sample_id)

    def __contains__(self, item):
        return item in self.__samples

    def __iter__(self):
        #iterate over samples of the group
        for sample in self.__samples:
            yield sample

    def __len__(self):
        #returns the number of samples of the group
        return len(self.__samples)

    def __str__(self):
        return str(self.__samples)

    def __repr__(self):
        return str(self)

class samples(object):
    """samples class describes the overall set of samples. Samples in different conditions
    are in different groups """
    #describes the set of samples. there are different samples in different conditions

    def __init__(self, list_samples, list_covariates, covariate_order=None):
        self.__groups = dict()
        self.__association = dict()
        self.__group_order = covariate_order
        # self.__load_samples(sample_filename)
        for sample, group in zip(list_samples, list_covariates):
            self.__add_sample_to_group(sample, group)
            self.__association[sample] = group

    def get_group(self, sample_id):
        return self.__association[sample_id]


    @property
    def samples(self):
        """It returns a list containing the samples' names.
        Samples of the same group are 'messi vicini' ;D  """
        #returns a list containing the sample names. Samples are grouped by group (xd)
        ret = list()
        #sort samples of the same group in lexicographic order
        grouped_samples = [
            (sorted(group.samples), group.group_id) for group in self.__groups.values()
        ]
        if self.__group_order:
            #sort groups using order given by the user
            grouped_samples = {group_id: group_list for group_list, group_id in grouped_samples}
            for curr_group in self.__group_order:
                ret.extend(grouped_samples[curr_group])
        else:
            #sort groups in lexicographic order
            for sublist, _ in sorted(grouped_samples, key=lambda t: t[1]):
                ret.extend(sublist)

        return ret

    def __add_sample_to_group(self, sample_id, group_id):
        if group_id not in self.__groups:
            self.__groups[group_id] = sample_group(group_id)

        self.__groups[group_id].add_sample(sample_id)

    def __contains__(self, item):
        #override in operator
        return any([item in sample_list for sample_list in self.__groups.values()])

    def __iter__(self):
        """Iterate over groups of samples. Returns sample_group objects"""
        for group in self.__groups.keys():
            yield self.__groups[group]


class circRNA(object):
    """Class circRNA describes a single circRNA.
    CircRNAs are identified by their coordinates: chromosome, start, end, strand. """

    #describes a single circRNA
    def __init__(self, circ_id):
        self.__id = circ_id         #circ_id: chr_start_end_strand
        self.__samples = dict()     #list of samples where the circRNA appears

    def add_sample(self, sample_id, num_reads):
        self.__samples[sample_id] = num_reads

    def get_reads_by_sample(self, sample_id):
        return self.__samples[sample_id] \
                if sample_id in self.__samples \
                else 0

    def check_min_replicates(self, min_reads, min_replicates, group_samples):
        #check if the circ has at least <min_replicates> replicates with
        #at least <min_reads> reads in at least one group.
        check = False

        for group in group_samples:
            nrep = 0
#            print("checking group {}".format(group.group_id))

            for sample in group:
                if sample in self.__samples and self.__samples[sample] >= min_reads:
#                    print("sample {} --> {}".format(sample, self.__samples[sample]))
                    nrep += 1

            if nrep >= min_replicates:
                check = True
                break

        return check

    def check_avg_read_number(self, threshold, group_samples):
        #check if the average number of reads in at least  group of samples is
        #greater of a given threshold
        check = False

        for group in group_samples:
            reads_per_sample = [self.__samples[sample] for sample in group if sample in self.__samples]
            avg_value = sum(reads_per_sample) / len(group)

            #print("group {} --> avg: {}".format(group.group_id, avg_value))
            if avg_value >= threshold:
                check = True
                break

        return check


    @property
    def circ_id(self):
        return self.__id

    @property
    def samples(self):
        return dict(self.__samples)

    def __str__(self):
        return self.__id + " " + str(self.__samples)



class set_circRNA(object):
    #describes a set of circRNA that appear in a set of samples
    def __init__(self, samples):
        self.__circs = dict()
        #self.__samples = set()
        self.__samples = samples    #type: <class 'samples'>


    def init_samples(self, sample_filename):
        #associate each sample with a group (using the provided file)
        with open(sample_filename) as fi:
            ficsv = csv.reader(fi, delimiter="\t")
            #skip header
            next(ficsv)

            for line in ficsv:
                samplename, group = line[:2]
                self.__samples.add_sample_to_group(samplename, group)

    def __iter__(self):
        for circ in self.__circs.values():
            yield circ

    def __len__(self):
        return len(self.__circs)

    def apply_replicate_filter(self, min_reads, min_replicates):
        #remove those circs that does not have at least <min_reads> reads
        #in <min_replicates> replicates in at least 1 group

        #lists circ ids of circs that does not pass the check
        deleted_circs = [circ.circ_id for circ in self \
            if not circ.check_min_replicates(min_reads, min_replicates, self.__samples)
        ]

#        print("Debug: ", len(deleted_circs))
        for circ_id in deleted_circs:
            del self.__circs[circ_id]


    def apply_average_filter(self, average):
        #remove those circs that does not have at least <average> reads in at least one group

        deleted_circs = [circ.circ_id for circ in self \
            if not circ.check_avg_read_number(average, self.__samples)
        ]

        for circ_id in deleted_circs:
            del self.__circs[circ_id]


    def get_circ(self, key="chr9_33948372_33948585_-"):
        return self.__circs[key]



    def add_circRNA(self, circ_id, sample, num_reads):
        if circ_id not in self.__circs:
            self.__circs[circ_id] = circRNA(circ_id)

        self.__circs[circ_id].add_sample(sample, num_reads)

    def write_circRNAs(self, filename, used_tool):
        """ This method creates two output files """

        sorted_samples = self.__samples.samples

        samples = {sample: 0 for sample in sorted_samples}

        with open("{}.{}.txt".format(filename, used_tool), "w") as crna_file, \
             open("{}.{}.count_table".format(filename, used_tool), "w") as crna_counts:
            #fowriter = csv.writer(fo, delimiter="\t")
            crnawriter = csv.writer(crna_file, delimiter="\t")
            crna_countswriter = csv.writer(crna_counts, delimiter="\t")

            #write header: column names are in the form SAMPLE_COVARIATE
            samples_covs = ["{}_{}".format(sample, self.__samples.get_group(sample)) for sample in sorted_samples]
            header = ["#circ_id"] + samples_covs
            crna_countswriter.writerow(header)

            #sort circs by id (i'd like to sort them by chromosome...)
            circs = [(id_circ, vals) for id_circ, vals in self.__circs.items()]
            circs.sort(key=lambda x: x[0])

            for circ_id, circ_object in circs:
                #write count file (statistical purposes)
                current = [circ_id]

                for sample_id in sorted_samples:
                    num_reads = circ_object.get_reads_by_sample(sample_id)
                    current.append(num_reads)

                crna_countswriter.writerow(current)

                #write circular rna file (CircHunter input file)
                chr, start, end, strand = circ_id.split("_")
                strand = 1 if strand == "+" else -1
                circ_name = "{}_{}_{}".format(chr, start, end)
                crnawriter.writerow([chr, start, end, circ_name, strand])



class CircRNA_Arguments(object):
    def __init__(self):
        #(chrm, start, end, strand, #reads)
        self.parameters = {
            "ciri":     
            (1, 2, 3, False, 4), 
            "ciri2":    
            (1, 2, 3, 10, 4),
            "acfs":     
            (0, 1, 2, 5, 4), 
            "circexplorer":
            (0, 1, 2, 5, 12), 
            "circexplorer2":
            (0, 1, 2, 5, 12), 
            "circrnafinder":
            (0, 1, 2, 5, 4), 
            "dcc":
            (0, 1, 2, 5, 3), 
            "findcirc2":
            (0, 1, 2, 5, 4), 
            "knife":
            (0, False, False, False, 1)
        }
    
    def __getitem__(self, key):
        if key in self.parameters:
            return self.parameters[key]
        return None 




def process_tool_line(tool, fields, line):
    """ Return a pair (circRNA_id, num_reads) where num_reads is the expression level of that circRNA """ 
    chrm, start, end, strand, num_reads = [line[i] if i else None for i in fields]

    if tool == "knife":
        chrm, start, end, strand = process_knife_line(line[chrm])
    elif tool == "ciri":
        strand = "+"
        
    chrm = chrm.lower()
    if "chr" not in chrm:
        chrm = "chr{}".format(chrm)

    circ_id = "{}_{}_{}_{}".format(chrm.replace("_", "-"), start, end, strand)
    return circ_id, int(num_reads)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #I/O
    parser.add_argument("-i", "--in", dest="input_dir", action="store", type=str, required=True)
    parser.add_argument("-o", "--out", dest="output_file", action="store", type=str, required=True)

    parser.add_argument("-s", "--samples", dest="samples", action="store", nargs="+", type=str, required=True)
    parser.add_argument("--cov", dest="covariates", action="store", nargs="+", type=str, required=True)
    parser.add_argument("--order", dest="order_by", action="store", nargs="+", type=str)
    parser.add_argument("--file", dest="file_extension", action="store", type=str, choices=list(CircRNA_Arguments().parameters.keys()), required=True)
    #filter parameters
    parser.add_argument("--mr", dest="min_reads", action="store", type=int, default=0)
    parser.add_argument("--mrep", dest="min_replicates", action="store", type=int, default=0)
    parser.add_argument("--avg", dest="average", action="store", type=float, default=10)

    args = parser.parse_args()

    if len(args.samples) != len(args.covariates):
        raise Exception("Error: samples and covariates lists must have the same length!")    

    useful_indexes = CircRNA_Arguments()[args.file_extension]

    samples_obj = samples(args.samples, args.covariates, args.order_by)
    circRNAs = set_circRNA(samples_obj)

    print("Reading sample files...")

    for sample_dir in os.listdir(args.input_dir):
        #we assume that each sample directory is named as the sample 
        if sample_dir in args.samples: 
            prefix = os.path.join(args.input_dir, sample_dir)
            #obtaining full filename of the sample 
            circRNA_prediction_file = [
                os.path.join(prefix, f) \
                for f in os.listdir(prefix) \
                if f.split(".")[-1].lower() == args.file_extension]#.pop()

            if len(circRNA_prediction_file) == 0:
                print("Skipping {} sample because its extension does not match with '.{}'".format(sample_dir, args.file_extension))
                continue

            circRNA_prediction_file = circRNA_prediction_file.pop()
            print("Processing {}".format(circRNA_prediction_file))

            with open(circRNA_prediction_file) as fi: 
                csvreader = csv.reader(fi, delimiter="\t")

                #test header presence 
                first_line = next(csvreader, None)
                if first_line is not None: 
                    #try to extract #num_reads and cast it to int
                    try: 
                        nreads = int(first_line[useful_indexes[-1]])
                        fi.seek(0, 0) 
                    except ValueError:
                        pass  #header spotted 

                for line in csvreader:
                    circ_id, circ_exp = process_tool_line(args.file_extension, useful_indexes, line)
                    circRNAs.add_circRNA(circ_id, sample_dir, circ_exp)


    print("{} circRNAs collected".format(len(circRNAs)))

    print("Applying filters...")

    circRNAs.apply_replicate_filter(args.min_reads, args.min_replicates)

    print("{} circRNAs passed the first check".format(len(circRNAs)))

    circRNAs.apply_average_filter(args.average)

    print("{} circRNAs passed the second check".format(len(circRNAs)))

    print("Writing output file...")
    circRNAs.write_circRNAs(args.output_file, args.file_extension)
