#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import re
import kkpdb
import kkstruct
from prepare_chain_terminus import detect_type_of_chain 
from prepare_chain_terminus import modify_canonical_resname
from optparse import OptionParser
from kkdefine import KKDEF as kk
import kkmdconf 

def define_options():
    p = OptionParser()
    p.add_option('-i', dest='fn_pdb',
                 help="")
    p.add_option('--i-orig', dest='fn_orig_pdb',
                 help="")
    p.add_option('--i-mdconf', dest='fn_mdconf',
                 help="")
    p.add_option('-o', dest='fn_out',
                 help="Outpu pdb file name")
    p.add_option('--o-pdb-canoresid', dest='fn_pdb_canoresid',
                 help="Output pdb file name with canonical residue id")
    p.add_option('--o-atom', dest='fn_atom',
                 help="")
    p.add_option('--o-res', dest='fn_res',
                 help="")
    p.add_option('--o-atom-woh', dest='fn_atom_woh',
                 help="without header line")
    p.add_option('--o-res-woh', dest='fn_res_woh',
                 help="without header line")
    p.add_option('--o-chain', dest='fn_chain',
                 help="")
    p.add_option('-s', dest="skip_ids",
                 type="int",
                 action="append",
                 default=[],
                 help="residues does not exist in original pdb")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = define_options() 
    if not opts.fn_pdb:
        return
    print "READ PDB: " + opts.fn_pdb
    reader = kkpdb.PDBReader(opts.fn_pdb)

    seg_undefined = "undefined"

    model = reader.read_model()
    model.set_residues_from_atom_info()
        
    mdconf = None
    if opts.fn_mdconf:
        mdconf = kkmdconf.MDConfReader(opts.fn_mdconf).read_conf()

    ##chains = model.split_with_chain()
    #for chain_id, chain in chains.items():
    #    modify_canonical_resname(chain)
    model.set_chain_types()
    #for chain_id, chain in model.chains.items():
    #    print chain_id + " : " + " ".join([str(x) for x in chain.res_indice])
    ## model.checking_disconnect_of_residues(reassign_chain_id=True)
    ##detect chain type
    #chain_types = {}
    #for chain_id in sorted(chains.key()):
    #    chain = chains[chain_id]
    #    chain_types[chain_id] = detect_type_of_chain(chain)
        ##model.set_residues_from_atom_info()

    reader_orig = None
    if opts.fn_orig_pdb:
        reader_orig = kkpdb.PDBReader(opts.fn_orig_pdb)
        model_orig = reader_orig.read_model()
        model.set_res_num_auth_from_model(model_orig, skip_ids=set(opts.skip_ids))
    model.set_chain_types()

    ##output atom info
    f_atom = open(opts.fn_atom,"w")
    f_atom_woh = None
    if opts.fn_atom_woh:
        f_atom_woh = open(opts.fn_atom_woh,"w")
    line = "\t".join(["atom_id.int","atom_name.string","res_name.string",
                      "res_id.int", "res_num.int", "res_num_auth.int",
                      "chain_id.string", "segment.string", "seg_type.string"])
    f_atom.write(line+"\n")
    for i,atom in enumerate(model.atoms):
        seg = seg_undefined
        seg_type = seg_undefined
        if mdconf:
            seg, seg_info = mdconf.get_segment_bynum(i, default=seg_undefined)
            seg_type = seg_info[2]

        line = "\t".join([str(x) for x in (atom.atom_id, atom.atom_name, atom.res_name, atom.res_index, atom.res_num, atom.res_num_auth, atom.chain_id, seg, seg_type)])
        f_atom.write(line+"\n")
        if f_atom_woh: f_atom_woh.write(line+"\n")
    f_atom.close()
    if f_atom_woh: f_atom_woh.close()
    ##output res info

    f_res= open(opts.fn_res,"w")
    f_res_woh = None
    if opts.fn_res_woh: f_res_woh = open(opts.fn_res_woh, "w")
    line = "\t".join(["res_id.int","res_num.int","res_num_auth.int","res_name.string",
                      "first_atom_id.int", "last_atom_id.int", "res_label.string",
                      "chain_id.string", "segment.string", "seg_type.string"])
    f_res.write(line+"\n")

    for res_index, residue in model.residues.items():
        res_label = residue.res_name[0].upper() + residue.res_name[1:].lower() + str(residue.res_num_auth)
        m = re.compile("d([atgc])", re.IGNORECASE).match(residue.res_name)
        if m:
            res_label = residue.res_name[1].upper() + str(residue.res_num_auth)
        seg = seg_undefined
        seg_type = seg_undefined
        if mdconf:
            seg, seg_info = mdconf.get_segment_bynum(min(residue.atom_indice), default=seg_undefined)
            seg_type = seg_info[2]
        chain_id = model.atoms[min(residue.atom_indice)].chain_id
        line = "\t".join([str(x) for x in (residue.res_index, residue.res_num, residue.res_num_auth, residue.res_name, min(residue.atom_indice), max(residue.atom_indice), res_label, chain_id, seg, seg_type)])
        f_res.write(line+"\n")
        if f_res_woh: f_res_woh.write(line+"\n")
    f_res.close()
    if f_res_woh: f_res_woh.close()

    ##output chain info
    f_chain = open(opts.fn_chain,"w")
    line = "\t".join(["chain_id.string","chain_type.string",
                      "first_atom_id.int", "last_atom_id.int",
                      "first_res_id.int", "last_res_id.int"
                      ])
    f_chain.write(line+"\n")
    for chain_id, chain in model.chains.items():
        #print chain_id
        #print chain.res_indice
        line = "\t".join([str(x) for x in (chain.chain_id, chain.chain_type, min(chain.atom_indice), max(chain.atom_indice),min(chain.res_indice), max(chain.res_indice))])
        f_chain.write(line+"\n")
    f_chain.close()

    writer = kkpdb.PDBWriter(opts.fn_out)
    if opts.fn_orig_pdb: model.swap_res_num_auth()
    writer.write_model(model)

    if opts.fn_pdb_canoresid:
        writer = kkpdb.PDBWriter(opts.fn_pdb_canoresid)
        model.renumber_res_num()
        writer.write_model(model)

    return



if __name__ == "__main__":
    _main()
