#!/usr/bin/python3

import argparse
import csv
import os, sys

# this script merges 2+ CIRI2 output files.
# it produces an output file containing for each circRNA,
# the number of junction reads of each sample
# Parameters:
#   1. input_folder that contains CIRI2 output files to merge
#   2. samplefile that contains AH NON LO SO IO
#   3. output name that will be the output filename



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

class samples(object):
    """Class samples describes the overall set of samples. Samples in different conditions
    are in different groups """
    #describes the set of samples. there are different samples in different conditions
    def __init__(self, sample_filename):
        self.__groups = dict()
        self.__load_samples(sample_filename)

    @property
    def samples(self):
        """It returns a list containing the samples' names.
        Samples of the same group are 'messi vicini' ;D  """
        #returns a list containing the sample names. Samples are grouped by group (xd)
        ret = list()
        for group in self.__groups.values():
            ret.extend(sorted(group.samples))

        return ret

    def __load_samples(self, sample_filename):
        """It loads samples from the input file and
        associate them to the corresponding group"""
        with open(sample_filename) as fi:
            ficsv = csv.reader(fi, delimiter="\t")
            #skip header
            next(ficsv)

            for line in ficsv:
                samplename, group = line[:2]
                self.__add_sample_to_group(samplename, int(group))


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

#            print("group {} --> avg: {}".format(group.group_id, avg_value))
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

    def write_circRNAs(self, filename):
        """ This method creates two output files """

        sorted_samples = self.__samples.samples
        samples = {sample: 0 for sample in sorted_samples}

        with open("{}.crna".format(filename), "w") as crna_file, \
             open("{}.crna_count".format(filename), "w") as crna_counts:
            #fowriter = csv.writer(fo, delimiter="\t")
            crnawriter = csv.writer(crna_file, delimiter="\t")
            crna_countswriter = csv.writer(crna_counts, delimiter="\t")

            #write header
            header = ["#circ_id"] + sorted_samples
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    #input/output parameters
    parser.add_argument("-i", "--in", dest="input_dir", action="store", type=str, required=True)
    parser.add_argument("-o", "--out", dest="output_file", action="store", type=str, required=True)
    parser.add_argument("-g", "--groups", dest="sample_group_file", action="store", type=str, required=True)

    #filter parameters
    parser.add_argument("--mr", dest="min_reads", action="store", type=int, default=2)
    parser.add_argument("--mrep", dest="min_replicates", action="store", type=int, default=2)
    parser.add_argument("--avg", dest="average", action="store", type=float, default=10)

    args = parser.parse_args()

    path = args.input_dir  #sys.argv[1]  #--in
    output = args.output_file #sys.argv[2]    #--out
    samplegroup_file = args.sample_group_file #sys.argv[3]  #--groups

    samples_obj = samples(samplegroup_file)
    circRNAs = set_circRNA(samples_obj)


    print("Reading sample files...", end="")

    for filename in os.listdir(path):
        filepath = os.path.join(path, filename)
        filename_noext = ".".join(filename.split(".")[:-1])

        if os.path.isfile(filepath) and filename.split(".")[-1].lower() == "ciri":
            sample = filename_noext.split("_")[0]
#            print("Processing sample {}... ".format(sample), end="")

            if sample in samples_obj:
                nlines = 0

                with open(filepath) as fi:
                    fireader = csv.reader(fi, delimiter="\t")
                    next(fireader)  #skip header

                    for line in fireader:
                        nlines += 1
                        chr = line[1] if "chr" in line[1].lower() else "chr{}".format(line[1])
                        circ_id = "{}_{}_{}_{}".format(chr, line[2], line[3], line[10])
                        num_reads = int(line[4])

                        circRNAs.add_circRNA(circ_id, sample, num_reads)

#                print("{} circRNAs processed".format(nlines))
#            else:
#                print("skipped".format(sample))

    print("{} circRNAs collected".format(len(circRNAs)))

    print("Applying filters...")

    circRNAs.apply_replicate_filter(args.min_reads, args.min_replicates)

    print("{} circRNAs passed the first check".format(len(circRNAs)))

    circRNAs.apply_average_filter(args.average)

    print("{} circRNAs passed the second check".format(len(circRNAs)))


    print("Writing output file...")
    circRNAs.write_circRNAs(output)
