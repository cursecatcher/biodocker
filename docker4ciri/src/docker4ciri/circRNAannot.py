#!/usr/bin/python3
# -*- coding: utf-8 -*-

from abc import ABC, abstractmethod      #abstract base classes
import argparse
import csv
import enum
import os
import pyliftover
import re
import tarfile
import urllib.request


class AssemblyVersion(enum.Enum):
    #human
    HG18 = "hg18"
    HG19 = "hg19"
    HG38 = "hg38"
    #mouse 
    MM9 = "mm9"
    MM10 = "mm10"

    @classmethod
    def has_value(cls, val):
        return any(val == item.value for item in cls)


class Organism(enum.Enum):
    "Supported organisms and their associated available assemblies "
    HUMAN = (AssemblyVersion.HG18, AssemblyVersion.HG19, AssemblyVersion.HG38)
    MOUSE = (AssemblyVersion.MM9, AssemblyVersion.MM10)

    @classmethod
    def get_organism_by_assembly(cls, v_assembly):
        """Return the organism associated to the given assembly version, if it is supported"""

        for organism in Organism: 
            if v_assembly in set(assembly.value for assembly in organism.value):
                return organism 
        raise Exception("Unrecognized assembly.")

    def __str__(self):
        return {Organism.HUMAN: "human", Organism.MOUSE: "mouse"}[self]



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
    """ Abstract class implementing a generic circRNA database """ 

    def __init__(self, available_dbs, assembly):
        self.__assembly = assembly
        self.__urls = available_dbs 

        is_supported = self.is_assembly_supported(assembly)
        if is_supported is None: 
            raise Exception("The selected db doesn't provide data about the given organism.")
        assembly_db, url = is_supported

        self.__url =  [url] if type(url) is str else url 
        self.__lifter = Lifter(assembly, assembly_db)

        self.__header = None
        self.__path_db = list() 
        self._annotations = list()


    def is_assembly_supported(self, my_assembly):
        """ Check if the organism associated to the given assembly has a db available """ 
        best_candidate = None

        if self.__urls is not None: 
            my_organism = Organism.get_organism_by_assembly(my_assembly)

            for assembly, url in self.__urls.items():
                current_organism = Organism.get_organism_by_assembly(assembly)

                if my_assembly == assembly:
                    best_candidate = assembly, url 
                    break 
                elif current_organism == my_organism:
                    best_candidate = assembly, url 

        return best_candidate  

    @abstractmethod
    def annotate(self, circset):
        pass

    def download_db(self, dest):
        """ Download the database(s) from URL and save the file in the specified directory """

        for db_url in self.__url: 
            print("downloading {}...".format(db_url))
            db_filename = os.path.join(dest, os.path.basename(db_url))
            #check if file has been already downloaded. Otherwise, it is downloaded now
            if not os.path.isfile(db_filename):
                print("Downloading db file from {}...".format(db_filename))
                response = urllib.request.urlopen(db_url)
                #downloaded_data = response.read() 
                with open(db_filename, "wb") as fdest:
                    fdest.write(response.read())
            #if downloaded file is an archive, it is extracted
            if ".tar.gz" in db_filename: #TODO -> db_filename.endswith(".tar.gz") ??
                with tarfile.open(db_filename, "r:gz") as tar: 
                    content = tar.getmembers()[0].name 
                    db_filename = os.path.join(dest, content) 

                    if not os.path.isfile(db_filename):
                        tar.extract(member=content, path=dest)
            
            self.__path_db.append(db_filename)


    def save(self, output_folder, output_file):
        """Save the annotated circRNAs in the specified destination """ 

        #create directory, if it does not exist
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        #obtain full path
        output_path = os.path.join(output_folder, output_file)
        #write annotation data
        with open(output_path, "w") as f:
            csvreader = csv.writer(f, delimiter="\t")

            if len(self._annotations) > 0:
                #check if header is present
                if self._annotations[0][0] is None:
                    csvreader.writerow(self._annotations[0][1])
                    self._annotations.pop(0)

                for circ, annotation in self._annotations:
                    csvreader.writerow([circ.id] + annotation)


    def __len__(self):
        """ Return the number of annotated circRNAs """
        return len(self._annotations)

    @property
    def path_db(self):
        return self.__path_db
    
    @path_db.setter
    def path_db(self, value):
        self.__path_db = value 

    @property
    def assembly(self):
        return self.__assembly
    
    @property
    def lifter(self):
        return self.__lifter


class CircBaseDB(circrnaDB):
    def __init__(self, assembly):
        urls = {
            AssemblyVersion.HG19.value: "http://www.circbase.org/download/hsa_hg19_circRNA.txt", 
            AssemblyVersion.MM9.value: "http://www.circbase.org/download/mmu_mm9_circRNA.bed"
        }
        super().__init__(urls, assembly)

    def annotate(self, circset):
        with open(self.path_db[0]) as f:
            annotation_file = csv.reader(f, delimiter="\t")

            indexes = None #chromosome, start, end and strand indexes in a line 

            if f.name.endswith(".txt"):
                #human version 
                header = next(annotation_file)    #skip header
                self._annotations.append((None, ["id"] + header[4:]))    #header
                indexes = 0, 1, 2, 3
            else:
                #mouse version 
                header = "TODO" #TODO - add some header plz
                indexes = 0, 1, 2, 4 

                # for line in annotation_file:
                #     chromo, start, end = line[:3]
                #     chromo, start, end, strand, flag = self.lifter.convert_coordinate(chromo, start, end, line[4])
                #     start += 1

                #     if flag and circset.check_circ(chromo, start, end, strand):
                #         circ = circRNA(chromo, start, end, strand)
                #         self._annotations.append((circ, line[4:]))
            

            for line in annotation_file:
                chromo, start, end, strand = [line[index] for index in indexes]
                chromo, start, end, strand, flag = self.lifter.convert_coordinates(chromo, start, end, strand)
                start += 1

                if flag and circset.check_circ(chromo, start, end, strand):
                    circ = circRNA(chromo, start, end, strand)
                    #obtain all the fields except coordinates 
                    fields = [line[index] for index in range(len(line)) if index not in indexes]
                    self._annotations.append((circ, fields))

        return self

    def download_db(self, output_folder):
        super(CircBaseDB, self).download_db(output_folder)
        return self 

    def save(self, output_folder):
        super(CircBaseDB, self).save(output_folder, "circbase.anno")
        return self 

class tscdDB(circrnaDB):
    def __init__(self, assembly):
        urls = {
            AssemblyVersion.HG19.value: (
                "http://gb.whu.edu.cn/TSCD/download/hg19_adult_TS_circRNAs.tar.gz", 
                "http://gb.whu.edu.cn/TSCD/download/hg19_fetal_TS_circRNAs.tar.gz"), 
            AssemblyVersion.HG38.value: (
                "http://gb.whu.edu.cn/TSCD/download/hg38_adult_TS_circRNAs.tar.gz", 
                "http://gb.whu.edu.cn/TSCD/download/hg38_fetal_TS_circRNAs.tar.gz")
        }
        
        super().__init__(urls, assembly)


#        super(tscdDB, self).save(output_folder, "tscd_{}.anno".format(self.__version))
    def annotate(self, circset):
        header = ["ID", "Sample_ID", "TSCD_version", "Junction", "Algorithm",
                    "BS_read", "Symbol", "CircRNA_type", "Region",
                    "Strand", "MRE", "RBP", "Genomic"]
        self._annotations.append((None, header))

        for path_db in self.path_db:
            with open(path_db, encoding="latin-1") as f:
                tscd_version = "fetal" if "fetal" in path_db else "adult"

                for line in csv.reader(f, delimiter="\t"):
                    chromo, start, end = line[2:5]
                    strand = line[6]

                    chromo, new_start, new_end, strand, flag = self.lifter.convert_coordinates(chromo, start, end, strand)
                    new_start += 1  #stupid different conventions

                    if flag and circset.check_circ(chromo, new_start, new_end, strand):
                        circ = circRNA(chromo, start, end, strand)
                        var = [line[0], line[5], tscd_version] + line[7:] 
                        self._annotations.append((circ, var))

        return self

    def download_db(self, output_folder):
        super(tscdDB, self).download_db(output_folder)
        return self 

    def save(self, output_folder):
        super(tscdDB, self).save(output_folder, "tscd.anno")
        return self 


class ExoRBaseDB(circrnaDB): 
    split_circrna_id_regex = re.compile(r"[:-]")

    def __init__(self, assembly):
        urls = {
            AssemblyVersion.HG38.value: "http://www.exorbase.org/exoRBase/download/download?file=Samples_combined_circRNA_RPM.txt"
        }
        super().__init__(urls, assembly)

 
    def annotate(self, circset):
        with open(self.path_db[0]) as f: 
            annotation_file = csv.reader(f, delimiter="\t")
            header = next(annotation_file)

            self._annotations.append((None, ["id", header[0]] + header[2:])) #TODO - 
            
            for line in annotation_file: 
                #example: chrY:13854443-13871798:+
                tokens = ExoRBaseDB.split_circrna_id_regex.split(line[1])
                chr, start, end = tokens[:3]
                strand = line[1][-1]
                other_fields = [line[0]] + line[2:]

                chr, start, end, strand, flag = self.lifter.convert_coordinates(chr, start, end, strand)
                start += 1 

                if flag and circset.check_circ(chr, start, end, strand):
                    circ = circRNA(chr, start, end, strand)
                    self._annotations.append((circ, other_fields))              
        
        return self
            
    def download_db(self, output_folder):
        super(ExoRBaseDB, self).download_db(output_folder)
        return self 

    def save(self, output_folder):
        super(ExoRBaseDB, self).save(output_folder, "exorbase.anno")
        return self 


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
