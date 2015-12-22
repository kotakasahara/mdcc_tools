
class Element:
    def __init__(self, name, atom_num, mass, radius):
        self.name = name
        self.atom_num = atom_num
        self.mass = mass
        self.radius = radius

class KKDEF(object):
##converting three-letter residue name into one-letter
    AA_3_1 = {"ALA":"A",          "ASP":"D",
              "ARG":"R",          "ASN":"N",
              "CYS":"C",          "GLU":"E",
              "GLN":"Q",          "GLY":"G",
              "HIS":"H",          "ILE":"I",
              "HISE":"H",         "HIS":"H",
              "LEU":"L",          "LYS":"K",
              "MET":"M",          "PHE":"F",
              "PRO":"P",          "SER":"S",
              "THR":"T",          "TRP":"W",
              "TYR":"Y",          "VAL":"V",
              "ACE":"*",          "NME":"*",
              "H1E":"H",          "H2E":"H",
              "H1D":"H",          "H2D":"H",
              "S1P":"S",          "S2P":"S",
              "T1P":"T",          "T2P":"T",
              "Y1P":"Y",          "Y2P":"Y",
              "DUM":"*"}
    NA_3_1 = {"GUA":"G",          "DG":"G",
              "DGU":"G",          "CYT":"C",
              "DC":"C",           "DCY":"C",
              "ADE":"A",          "DA":"A",
              "DAD":"A",          "THY":"T",
              "DT":"T",           "DTH":"T",
              "DUM":"*"
              }
    
    RESNAMES = {"ALA":"ALA",          "ASP":"ASP",
                "ARG":"ARG",          "ASN":"ASN",
                "CYS":"CYS",          "GLU":"GLU",
                "GLN":"GLN",          "GLY":"GLY",
                "HIS":"HIS",          "ILE":"ILE",
                "HISE":"HIS",         "HIE":"HIS",
                "LEU":"LEU",          "LYS":"LYS",
                "MET":"MET",          "PHE":"PHE",
                "PRO":"PRO",          "SER":"SER",
                "THR":"THR",          "TRP":"TRP",
                "TYR":"TYR",          "VAL":"VAL",

                ## terminal cap
                "ACE":"ACE",          "NME":"NME",

                ## phospho-amino acids
                "H1D":"H1D",          "H2D":"H2D",
                "H1E":"H1E",          "H2E":"H2E",
                "T1P":"T1P",          "T2P":"T2P",
                "S1P":"S1P",          "S2P":"S2P",
                "Y1P":"Y1P",          "Y2P":"Y2P",

                ## nucleotides 
                "GUA":"DG",           "DG":"DG",
                "DGU":"DG",           "CYT":"DC",
                "DC":"DC",            "DCY":"DC",
                "ADE":"DA",           "DA":"DA",
                "DAD":"DA",           "THY":"DT",
                "DT":"DT",            "DTH":"DT",

                ## waters
                "HOH":"HOH", "SOL":"HOH", "WAT":"HOH", # waters
                "TIP":"HOH",

                ## ions
                "CIP":"CIP", "CIM":"CIM",
                "NA":"NA",   "K":"K",  "CL":"CL",

                "DUM":"*"

              }
    
    ELEM = {1:Element("H", 1, 1.00794, 1.2),
            6:Element("C", 6, 12.0107, 1.7),
            7:Element("N", 7, 14.0067, 1.55),
            8:Element("O", 7, 15.9994, 1.52),
            15:Element("P", 15, 30.973762, 1.8),
            16:Element("S", 16, 32.065, 1.8),
    }
    ##<DEFINE_CONSTANT>
    ##Type of chains
    CH_DUMMY = "UNKNOWN"
    CH_PEPTIDE = "PEPTIDE"
    CH_DNA = "DNA"
    CH_RNA = "RNA"
    CH_CHEMICAL = "CHEMICAL"
    ##</DEFINE_CONSTANT>

class PrestoDef:
    Atom_type = {
        "C":1,
        "CC":1,
        "CW":1,
        "CR":1,
        "CA":1,
        "CB":1,
        "C*":1,
        "CD":1,
        "CK":1,
        "CM":1,
        "CN":1,
        "CQ":1,
        "CV":1,
        "CX":1,
        "CY":1,
        "CZ":1,
        "CT":2,
        "H":3,
        "HO":4,
        "HS":5,
        "HC":6,
        "H1":7,
        "H2":8,
        "H3":   9,
        "HP":  10,
        "HA":  11,
        "H4":  12,
        "H5":  13,
        "HW":  14,
        "N":  15,
        "NA":  15,
        "N*":  15,
        "N2":  15,
        "N3":  15,
        "NB":  15,
        "NC":  15,
        "NP":  15,
        "NO":  15,
        "NY":  15,
        "NT":  15,
        "O":  16,
        "O2":  17,
        "OW":  18,
        "OH":  19,
        "OS":  20,
        "P":  21,
        "S":  22,
        "SH":  23,
        "IM":  24,
        "IP":  25,
        "K":  26,
        "Li":  27,          
        "Na":  28,
        "RA":  29,
        "CS":  30,
        "MG":  31,
        "C0":  32,
        "Zn":  33,
        "F":  34,
        "CL":  35,
        "BR":  36,
        "I":  37,
        "IB":  38,
        "LP":  39,
        "SE":  40
    }



        
