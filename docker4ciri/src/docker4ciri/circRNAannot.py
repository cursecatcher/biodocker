#!/usr/bin/python3
# -*- coding: utf-8 -*-

from abc import ABC, abstractmethod      #abstract base classes
import argparse
import csv
import enum
import os
import pyliftover
import tarfile
import urllib.request


class AssemblyVersion(enum.Enum):
    HG18 = "hg18"
    HG19 = "hg19"
    HG38 = "hg38"

    @classmethod
    def has_value(cls, val):
        return any(val == item.value for item in cls)


class Lifter(pyliftover.LiftOver):
    def __init__(self, src, dest):
        self.__convert = src != dest
        if self.__convert:
            super().__init__(src, dest)

    def convert_coordinates(self, chromosome, start, end, strand):
        """ """
        start, end = int(start), int(end)
        flag = True

        if self.__convert:
            new_start = super().convert_coordinate(chromosome, start, strand)
            new_end = super().convert_coordinate(chromosome, end, strand)
            flag = False

            if new_start is not None and new_end is not None:
                if len(new_start) == len(new_end) == 1:
                    start = new_start[0][1]
                    end = new_end[0][1]
                    flag = True

        return chromosome, start, end, strand, flag


class circrnaDB(ABC):
    def __init__(self, url):
        self._annotations = list()
        self.__url = url
        self.__filename = None
        self.__header = None

    @abstractmethod
    def annotate(self, circset):
        pass

    def download_db(self, dest):
        """ """
        filename = os.path.basename(self.__url)
        self.__filename = os.path.join(dest, filename)
        #check if file is already present. Otherwise, it is downloaded
        if not os.path.isfile(self.__filename):
            print("Downloading file from {}...".format(self.__url))

            response = urllib.request.urlopen(self.__url)
            downloaded_data = response.read()

            with open(self.__filename, "wb") as fdest:
                fdest.write(downloaded_data)
        #if downloaded file is an archive, it is extracted
        if ".tar.gz" in filename:
            with tarfile.open(self.__filename, "r:gz") as tar:
                content = tar.getmembers()[0].name
                self.__filename = os.path.join(dest, content)

                if not os.path.isfile(self.__filename):
                    tar.extract(member = content, path = dest)

        return self.__filename

    def save(self, output_folder, output_file):
        #create directory, if it does not exist
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        #obtain full path
        output_path = os.path.join(output_folder, output_file)
        #write annotation data
        with open(output_path, "w") as f:
            csvreader = csv.writer(f, delimiter="\t")

            #check if header is present
            if self._annotations[0][0] is None:
                csvreader.writerow(self._annotations[0][1])
                self._annotations.pop(0)

            for circ, annotation in self._annotations:
                csvreader.writerow([circ.id] + annotation)

    def __len__(self):
        return len(self._annotations)


class CircBaseDB(circrnaDB):
    url_circbase = "http://www.circbase.org/download/hsa_hg19_circRNA.txt"

    def __init__(self, assembly):
        if not AssemblyVersion.has_value(assembly):
            raise Exception("UNSUPPORTED ASSEMBLY VERSION")

        super().__init__(CircBaseDB.url_circbase)
        self.__path_db = None
        self.__assembly = assembly

    def annotate(self, circset):
        lifter = Lifter(self.__assembly, "hg19")

        with open(self.__path_db) as f:
            annotation_file = csv.reader(f, delimiter="\t")

            header = next(annotation_file)    #skip header
            self._annotations.append((None, ["id"] + header[4:]))    #header

            for line in annotation_file:
                chromo, start, end, strand, flag = lifter.convert_coordinates(*line[:4])
                start = start + 1

                if flag and circset.check_circ(chromo, start, end, strand):
                    circ = circRNA(chromo, start, end, strand)
                    self._annotations.append((circ, line[4:]))

        return self

    def download_db(self, output_folder):
        self.__path_db = super(CircBaseDB, self).download_db(output_folder)

    def save(self, output_folder):
        super(CircBaseDB, self).save(output_folder, "circbase.anno")

class tscdDB(circrnaDB):
    hg19_adult = "http://gb.whu.edu.cn/TSCD/download/hg19_adult_TS_circRNAs.tar.gz"
    hg19_fetal = "http://gb.whu.edu.cn/TSCD/download/hg19_fetal_TS_circRNAs.tar.gz"
    hg38_adult = "http://gb.whu.edu.cn/TSCD/download/hg38_adult_TS_circRNAs.tar.gz"
    hg38_fetal = "http://gb.whu.edu.cn/TSCD/download/hg38_fetal_TS_circRNAs.tar.gz"

    def __init__(self, assembly, version):
        if not AssemblyVersion.has_value(assembly):
            raise Exception("UNSUPPORTED ASSEMBLY VERSION")

        version = version.lower()

        if version not in ["adult", "fetal"]:
            raise Exception("Invalid TSCD version")

        url = tscdDB.hg19_adult if version == "adult" else tscdDB.hg19_fetal
        self.__db_assembly = "hg19" #used to initialize LiftOver object

        if assembly == "hg38":
            url = tscdDB.hg38_adult if version == "adult" else tscdDB.hg38_fetal
            self.__db_assembly = "hg38"

        super().__init__(url)
        self.__version = version
        self.__path_db = None
        self.__assembly = assembly


    def annotate(self, circset):
        #convert from your assembly version to db assembly version
        lifter = Lifter(self.__assembly, self.__db_assembly)

        with open(self.__path_db, encoding="latin-1") as f:
            annotation_file = csv.reader(f, delimiter="\t")

            header = ["ID", "Sample_ID", "Junction", "Algorithm",
                "BS_read", "Symbol", "CircRNA_type", "Region",
                "Strand", "MRE", "RBP", "Genomic"]
            self._annotations.append((None, header))


            for line in annotation_file:
                chromo, start, end = line[2:5]
                strand = line[6]

                chromo, new_start, new_end, strand, flag = lifter.convert_coordinates(chromo, start, end, strand)

                new_start += 1  #stupid different conventions

                if flag and circset.check_circ(chromo, new_start, new_end, strand):
                    circ = circRNA(chromo, start, end, strand)
                    var = [line[0], line[5]] + line[7:]
                    self._annotations.append((circ, var))

        return self

    def download_db(self, output_folder):
        self.__path_db = super(tscdDB, self).download_db(output_folder)

    def save(self, output_folder):
        super(tscdDB, self).save(output_folder, "tscd_{}.anno".format(self.__version))


class circRNA(object):
    def __init__(self, chromosome, start, end, strand):
        self.__id = "{}_{}_{}_{}".format(chromosome, start, end, strand)

        self.__chromosome = chromosome
        self.__start = start
        self.__end = end
        self.__strand = strand

    @property
    def id(self):
        return self.__id

    @property
    def chromosome(self):
        return self.__chromosome

    @property
    def start(self):
        return self.__start

    @property
    def end(self):
        return self.__end

    @property
    def strand(self):
        return self.__strand

    def __hash__(self):
        return hash(self.__id)

    def __str__(self):
        return self.__id


class CircRNAs(object):
    def __init__(self):
        self.__circs = set()

    def add_circ(self, chromosome, start, end, strand):
        self.__circs.add(circRNA(chromosome, start, end, strand).id)

    def check_circ(self, chromosome, start, end, strand):
        return circRNA(chromosome, start, end, strand).id in self.__circs
