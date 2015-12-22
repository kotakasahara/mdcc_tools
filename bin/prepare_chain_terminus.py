#!/usr/bin/python2.7

from optparse import OptionParser
import numpy
import sys
import kkpdb
import collections

def define_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_pdb',
                 help="file name for input PDB")
    p.add_option('-o', dest='fn_out',
                 help="file name for output PDB")
    p.add_option('-c', dest='chaintype',
                 action="append",
                 help="type of chain, D(na) or P(rotein). A:P, for chain A as protein.")
    p.add_option('--no-nme', dest='no_nme',
                 action="store_true",
                 help="capping the n-terminus by NME")
    p.add_option('--no-ace', dest='no_ace',
                 action="store_true",
                 help="capping the c-terminus by ACE")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

AA_RES = ["ALA","ASP","ARG","ASN","CYS",
          "GLU","GLN","GLY","HIS","ILE",
          "LEU","LYS","MET","PHE","PRO",
          "SER","THR","TRP","TYR","VAL",
          "ACE", "NME", "HIE",
          "H1D", "H2D", "H1E", "H2E",
          "S1P", "S2P", "T1P", "T2P",
          "Y1P", "Y2P", "CYM", "CYSZ",
          "HIEZ"]
NA_RES = ["DT","DA","DG","DC"]
WAT_RES = ["HOH","WAT","SOL"]

CH_PEPTIDE = 0
CH_DNA = 1
CH_OTHER = 2

def modify_canonical_resname(chain):
    for atom in chain.atoms:
        atom.res_name = atom.res_name.replace("GUA","DG")
        atom.res_name = atom.res_name.replace("THY","DT")
        atom.res_name = atom.res_name.replace("CYT","DC")
        atom.res_name = atom.res_name.replace("ADE","DA")
        if atom.atom_name == "C5A":
            atom.atom_name = "C7"
    return chain

def detect_type_of_chain(chain):
    n_aa = 0
    n_na = 0
    other = set()
    ch_type = CH_OTHER
    for atom in chain.atoms:
        if atom.res_name in AA_RES: n_aa+=1
        elif atom.res_name in NA_RES: n_na+=1
        else:
            print "UNKNOWN RESIDUE: " + atom.res_name
            other.add(atom.res_name)
    print "n_aa:" + str(n_aa)
    print "n_na:" + str(n_na)
    print "n_other:" + str(len(other))

    if n_aa > 0 and n_na == 0 and len(other) == 0:
        ch_type = CH_PEPTIDE
    elif n_aa == 0 and n_na > 0 and len(other) == 0:
        ch_type = CH_DNA
    else:
        print other
    return ch_type

def prepare_chain_peptide(chain, no_nme, no_ace):
    pop_atoms = []
    for i,atom in enumerate(chain.atoms):
        if not no_ace and atom.res_id == chain.atoms[0].res_id:
            if atom.atom_name == "C" or \
                    atom.atom_name == "O":
                atom.res_name = "ACE"
                print "ACE:" + atom.info_txt()
            elif atom.atom_name == "CA":
                atom.atom_name = "CH3"
                atom.res_name = "ACE"                
            else:
                pop_atoms.append(i)
                print "POPPED:" + atom.info_txt()    
        if not no_nme and atom.res_id == chain.atoms[-1].res_id:
            if atom.atom_name == "N":
                atom.res_name = "NME"                
                print "NME:" + atom.info_txt()
            elif atom.atom_name == "CA":
                atom.res_name = "NME"
                atom.atom_name = "CH3"
                print "NME:" + atom.info_txt()
            else:
                pop_atoms.append(i)
                print "POPPED:" + atom.info_txt()  
    chain.pop_atom(pop_atoms)
    return chain

def prepare_chain_dna(chain):
    pop_atoms = []
    for i,atom in enumerate(chain.atoms):
        if atom.res_id == chain.atoms[0].res_id:
            if atom.atom_name == "P" or \
                    atom.atom_name == "OP1" or \
                    atom.atom_name == "OP2":
                pop_atoms.append(i)
                print "POPPED:" + atom.info_txt()
        if atom.atom_name == "C5A":
            atom.atom_name = "C7"
    chain.pop_atom(pop_atoms)
    return chain

def prepare_chain(chain, no_nme, no_ace, ch_type=-1):
    chain = modify_canonical_resname(chain)
    if ch_type==-1: ch_type = detect_type_of_chain(chain)
    print "CHAIN " + chain.atoms[0].chain_id + " : TYPE " + str(ch_type)
    if ch_type == CH_PEPTIDE:
        chain = prepare_chain_peptide(chain, no_nme, no_ace)
    elif ch_type == CH_DNA:
        chain = prepare_chain_dna(chain)
    return chain

def _main():
    opts, args = define_options() 
    if not opts.fn_pdb:
        return

    print "ABOUT THIS PROGRAM:"
    print "For peptide:"
    print "The first residue changed into ACE"
    print "that means delete atoms other than CA, C, and O,"
    print "and chainging residue name to ACE."
    print "Peptide chains must be described in PDB by order of N to C"
    print ""
    print "For DNA:"
    print "If there are P, O1P, or O2P atoms exist, delete them."
    print "DNA chains must be described in PDB by order of 5' to 3'"
    print ""

    print "READ PDB: " + opts.fn_pdb
    reader = kkpdb.PDBReader(opts.fn_pdb)
    model = reader.read_model()
    #model.set_water_chain()
    #model.set_chain_types()

    model.checking_disconnect_of_residues(reassign_chain_id=False)
    
    chains = model.split_with_chain()
    print sorted(chains.keys())

    chaintypes = {}
    for v in opts.chaintype:
        tmp = v.split(":")
        if tmp[1] == "P": 
            chaintypes[tmp[0]] = CH_PEPTIDE
        elif tmp[1] == "D":
            chaintypes[tmp[0]] = CH_DNA
            
    prepared = []
    for chain_id in sorted(chains.keys()):
        chain = chains[chain_id]
        chtype = -1
        if chain_id in chaintypes: chtype=chaintypes[chain_id]
        prepared.append(prepare_chain(chain, opts.no_nme, opts.no_ace, chtype))
    for chain in prepared[1:]:
        prepared[0].push_atoms(chain.atoms)
    
    writer = kkpdb.PDBWriter(opts.fn_out)
    writer.write_model(prepared[0])

if __name__ == "__main__":
    _main()
