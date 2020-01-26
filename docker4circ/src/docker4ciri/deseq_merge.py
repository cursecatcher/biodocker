#!/usr/bin/python3

import argparse
import csv
import logging
import os, sys

#represent a single row of the merged file
class Record(object):
    """Represent a single line of the merged file.  """

    def __init__(self, record_id, num_samples):
        self.__id = record_id
        self.__samples_list = [0] * num_samples

    @property
    def list(self):
        return [self.__id] + self.__samples_list

    def __setitem__(self, key, value):
        self.__samples_list[key] = value

    def __getitem__(self, key):
        return self.__samples_list[key]

    def __hash__(self):
        return hash(self.__id)

    def __str__(self):
        return str(self.__samples_list)

class SamplesCollection(object):
    """Represent a collection of record from different samples. """

    def __init__(self, samples, covariates, covariate_order=None):
        #mapping from record_ids to records
        self.__records = dict()
        #mapping from samples to indexes. (sorted list of samples, grouped by covariate)
        self.__samples = dict()
        #mapping sample -> covariate
        self.__sample_cov_mapping = {sample: cov for sample, cov in zip(samples, covariates)}
        #mapping sample -> index in a record (give order)
        self.__indexes = dict()     #

        samples_per_cov = dict()

        #group samples respect to covariate and set the given order to covariates
        for sample, cov in zip(samples, covariates):
            #group samples
            if cov not in samples_per_cov:
                samples_per_cov[cov] = list()
            samples_per_cov[cov].append(sample)
        else:
            #set default order
            if covariate_order is None:
                covariate_order = sorted([x for x in set(covariates)])

            i = 0
            #sort groups of the same covariate
            for covariate in covariate_order:
                for sample in sorted(samples_per_cov[covariate]):
                    self.__samples[sample] = i
                    i += 1

    def get_header(self):
        sorted_samples = sorted(self.__samples.items(), key=lambda pair: pair[1])
        return ["{}_{}".format(sample, self.__sample_cov_mapping[sample]) for sample, _ in sorted_samples]


    def __len__(self):
        return len(self.__records)

    def add_record(self, sample, key, value):
        if key not in self.__records:
            self.__records[key] = Record(key, len(self.__samples))

        self.__records[key][self.__samples[sample]] = value



    def write(self, filename):
        with open(filename, "w") as count_file:
            writer = csv.writer(count_file, delimiter="\t")

            header = ["ID"] + self.get_header()
            writer.writerow(header)

            for record in self.__records.values():
                writer.writerow(record.list)




if __name__ == "__main__":
    logging.getLogger().setLevel(logging.WARNING)

    parser = argparse.ArgumentParser()

    #input/output parameters
    parser.add_argument("-i", "--in", dest="input_dir", action="store", type=str, required=True)
    parser.add_argument("-o", "--out", dest="output_file", action="store", type=str, required=True)

    parser.add_argument("-s", "--samples", dest="samples", action="store", nargs="+", type=str, required=True)
    parser.add_argument("--cov", dest="covariates", action="store", nargs="+", type=str, required=True)
    parser.add_argument("--order", dest="order_by", action="store", nargs="+", type=str)

    parser.add_argument("--col", dest="column", action="store", type=int, required=True)
    parser.add_argument("--ext", dest="extension", action="store", type=str, required=True)

    args = parser.parse_args()

    if len(args.samples) != len(args.covariates):
        raise Exception("Error: samples and covariates lists must have the same length!")

    path = args.input_dir
    output = args.output_file


    #init SamplesCollection
    collection = SamplesCollection(args.samples, args.covariates, args.order_by)

    for sample_dir in os.listdir(path):
        #assume samples directories are called as the samples
        if sample_dir in args.samples:
            prefix = os.path.join(path, sample_dir)
            #get absolute path of the file we are searching for.
            #we also assume that each directory contains one and only one interesting file
            pathfile = [os.path.join(prefix, f) for f in os.listdir(prefix) if f.split(".")[-1] == args.extension].pop()

            print("Processing {}".format(pathfile))

            with open(pathfile) as fi:
                reader = csv.reader(fi, delimiter="\t")
                next(reader)

                for line in reader:
                    row_id, row_val = line[0], line[args.column - 1]
                    collection.add_record(sample_dir, row_id, row_val)

    print("{} elements collected".format(len(collection)))

    print("Header: {}".format(collection.get_header()))

    print("Writing output file")
    collection.write(output)
