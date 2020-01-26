#!/usr/bin/env python3

import enum 
import pyliftover


### Some utilities 
## Liftover  
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


### Enums
class SupportedTool(enum.Enum):
    ACFS = "acfs"
    CIRI = "ciri"
    CIRI2 = "ciri2"
    STARCHIP = "starchip"
    CIRCEXPLORER = "circexplorer"
    CIRCEXPLORER2 = "circexplorer2"
    UROBORUS = "uroborus"
    CIRCRNAFINDER = "circrnafinder"
    FINDCIRC2 = "findcirc2"
    KNIFE = "knife"
    DCC = "dcc"
    

    @classmethod
    def get_enum(cls, str_value):
        if type(str_value) is str:
            for enum_value in SupportedTool:
                if enum_value.value == str_value:
                    return enum_value
        elif type(str_value) is SupportedTool:
            return str_value
        return None 


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
    
    @classmethod
    def get_enum_value(cls, str_value):
        if type(str_value) is str:
            for enum_value in AssemblyVersion:
                if enum_value.value == str_value:
                    return enum_value
        elif type(str_value) is AssemblyVersion:
            return str_value
        return None 


class Organism(enum.Enum):
    "Supported organisms and their associated available assemblies "
    HUMAN = (AssemblyVersion.HG18, AssemblyVersion.HG19, AssemblyVersion.HG38)
    MOUSE = (AssemblyVersion.MM9, AssemblyVersion.MM10)

    @classmethod
    def get_organism_by_assembly(cls, v_assembly):
        """Return the organism associated to the given assembly version, if it is supported"""

        organism = [o for o in Organism if v_assembly in o.value]

        if len(organism) == 0:
            raise UnsupportedOrganism(v_assembly)
        return organism.pop()

    def __str__(self):
        return {Organism.HUMAN: "human", Organism.MOUSE: "mouse"}[self]




## Exceptions 
class NoDatabaseException(Exception):
    def __init__(self, db_name, assembly):
        super().__init__("{} database doesn't provide annotation for {} assembly".format(db_name, assembly))

class UnsupportedOrganism(Exception):
    def __init__(self, assembly):
        super().__init__("Unrecognized assembly: {}".format(assembly))



def process_knife_line(id_field):
    """ Return a 4-tuple containing (chrm, start, end, strand) """ 
    # i.e. chr13|DNAJC3:96377506|DNAJC3:96375496|rev|+

    return tuple(
        token if ":" not in token else token.split(":")[1]
        for i, token in enumerate(id_field.split("|")) if i != 3
    )