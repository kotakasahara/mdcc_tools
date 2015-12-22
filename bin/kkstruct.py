#!/usr/bin/python2.7

import numpy
import copy
from kkdefine import KKDEF as kk
import sys
import re
DEBUG = True

class Atom(object):
    def __init__(self, header="ATOM  ", atom_id=-1, atom_name="DUM", res_name="DUM",
                 chain_id="X", res_num=-1, x=0.0, y=0.0, z=0.0,
                 ocp=1.0, tf=0.0, seg_id=0.0,
                 elem="X", chg=0.0, atom_id_pdb=-1, res_num_auth=-1, entity_id=-1):
        self.header = header
        self.atom_id = int(atom_id)
        self.atom_id_pdb = atom_id_pdb
        self.atom_name = atom_name
        self.res_name = res_name
        self.chain_id = chain_id
        self.chain_type = kk.CH_DUMMY
        self.res_num = res_num
        self.res_num_auth = res_num_auth
        self.res_index = -1
        self.crd = numpy.array([x,y,z])
        self.ocp = ocp
        self.tf = tf
        self.seg_id = seg_id
        self.elem = elem
        self.chg = chg
        self.alt_loc = " "
        self.entity_id = entity_id
        ## unique id of residue in model
        ## this value was originally defined in this program
        self.res_index = 0
        
    def distTo(self, other):
        return numpy.sqrt(numpy.sum(numpy.power((self.crd-other.crd),2)))
    def info_txt(self):
        info = ""
        info += str(self.header)  + "\t"
        info += str(self.atom_id)  + "\t"
        info += self.atom_name  + "\t"
        info += self.res_name  + "\t"
        info += self.chain_id  + "\t"
        info += str(self.res_num)  + "\t"
        info += str(self.crd[0])  + "\t"
        info += str(self.crd[1])  + "\t"
        info += str(self.crd[2])  + "\t"
        info += str(self.ocp)  + "\t"
        info += str(self.tf)   + "\t"
        info += self.seg_id  + "\t"
        info += self.elem  + "\t"
        info += str(self.chg)  + "\t"
        return info

class Entity(object):
    def __init__(self, entity_id, weight, description, ent_cate_type):
        self.entity_id = entity_id
        self.weight = weight
        description = re.sub(r'\s','_',description)
        self.description = description
        self.ent_cate_type = ent_cate_type
        self.ent_type = ent_cate_type
        self.ent_seq = "-"
        return 
    def set_type(self, ent_type):
        self.ent_type = ent_type
        return
    def set_seq(self, seq):
        self.ent_seq = re.sub(r'\n', '', seq)
        return

class AtomSubset(object):
    def __init__(self):
        self.atom_indice = set()
        return
    def push_atom_index(self,atom_index):
        self.atom_indice.add(atom_index)

class Chain(AtomSubset):
    def __init__(self, chain_id):
        super(Chain, self).__init__()
        self.chain_id = chain_id
        self.chain_type = kk.CH_DUMMY
        self.res_indice = set()
        self.entity_id = -1
    def push_res_index(self,res_index):
        self.res_indice.add(res_index)
        
class Residue(AtomSubset):
    def __init__(self, res_index, res_num, res_name, chain_id, res_num_auth=0):
        super(Residue, self).__init__()
        self.res_num = res_num
        self.res_index = res_index
        self.res_name = res_name
        self.res_num_auth = res_num_auth
        self.chain_id = chain_id

class Model(object):
    def __init__(self):
        self.atoms = []
        self.cryst = numpy.array([0.0,0.0,0.0])
        self.model_id = 0
        self.title = ""
        ## self.residues[residue_num] = Residue()
        self.residues = {}
        ## self.cahins[chain_id] = Chain()
        self.chains = {}
        self.title = ""
        self.pbc_box = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.pdb_id = "XXXX"
        self.entities = {}
    ### Subset ####
    def sub_model(self, atom_ids):
        new_m = Model()
        new_m.cryst = self.cryst
        for aid in atom_ids:
            new_m.push_atom(self.atoms[aid])
        return new_m

    ## subset info: res, chain
    def reset_chains_from_atom_info(self):
        self.chains = {}
        for i, atom in enumerate(self.atoms):
            chain_id = atom.chain_id
            if chain_id == " ": chain_id="X"
            if not chain_id in self.chains:
                self.chains[chain_id] = Chain(chain_id)
            self.chains[chain_id].push_res_index(atom.res_index)
            self.chains[chain_id].push_atom_index(atom.atom_id)
        return
    def set_residues_from_atom_info(self):
        new_res_index = -1
        pointer_res_name = ""
        pointer_res_num = 0
        new_chain_id = "@"
        pointer_chain_id = "-"
        pointer_chain_type = ""
        for i, atom in enumerate(self.atoms):
            if atom.chain_id != pointer_chain_id or \
                    atom.chain_type != pointer_chain_type:
                #new_chain_id = chr(ord(new_chain_id)+1)
                new_chain_id = atom.chain_id
                pointer_chain_type = atom.chain_type
                pointer_chain_id = atom.chain_id
                self.chains[new_chain_id] = Chain(new_chain_id)
                self.chains[new_chain_id].chain_type = pointer_chain_type
            if atom.res_num != pointer_res_num or \
                    atom.res_name != pointer_res_name or \
                    atom.chain_id != pointer_chain_id:
                new_res_index += 1
                self.residues[new_res_index] = Residue(new_res_index, atom.res_num, atom.res_name,
                                                       atom.chain_id, atom.res_num_auth)

                pointer_res_num = atom.res_num
                pointer_res_name = atom.res_name
            self.atoms[i].res_index = new_res_index
            self.residues[new_res_index].push_atom_index(i)
            self.chains[new_chain_id].push_atom_index(i)
            self.chains[new_chain_id].push_res_index(new_res_index)
        return self.atoms
    def set_chain_types(self, reassign_chain_id=False):
        for chain_id in sorted(self.chains.keys()):
            chain = self.chains[chain_id]
            types = set()
            for atom_index in chain.atom_indice:
                atom = self.atoms[atom_index]
                chain.entity_id = atom.entity_id
                tmp_type = kk.CH_DUMMY
                if atom.res_name[0:3] in kk.AA_3_1:
                    tmp_type = kk.CH_PEPTIDE
                elif atom.res_name in kk.NA_3_1:
                    tmp_type = kk.CH_DNA
                atom.chain_type = tmp_type
                types.add(tmp_type)
                ##print chain_id + " " + atom.res_name + " " + tmp_type
            if len(types) == 1:
                chain.chain_type = list(types)[0]
            else:
                chain.chain_type = kk.CH_DUMMY
            #print "set_chain_types CHAIN_ID:" + chain_id + " " + str(chain.chain_type)
    ##
    def split_with_chain(self, ignore_res=[]):
        """
        Generating subsets of this models for each chain
        Returning a dictionary of models.
        chain_models[chain_id] = Model
        """
        chain_atoms = {}
        for i,atom in enumerate(self.atoms):
            chid = atom.chain_id
            #if chid == " ": chid = "X"
            if atom.res_name in ignore_res or \
                    (atom.res_name in kk.RESNAMES and \
                         kk.RESNAMES[atom.res_name] in ignore_res):
                chid = " "
                atom.chain_id = " "
            if not atom.chain_id in chain_atoms:
                chain_atoms[chid] = []
            #print str(i) + " " + atom.atom_name + " [" + chid + "]"
            chain_atoms[chid].append(i)
        chain_models = {}
        for chain_id, atom_indices in chain_atoms.items():
            chain_models[chain_id] = self.sub_model(atom_indices)
        #for chain_id, chain in self.chains.items():
        #    chain_models[chain_id] = self.sub_model(chain.atom_indice)
        return chain_models
    ### manipurate Atoms ####

    def push_atom(self, atom):
        self.atoms.append(copy.deepcopy(atom))
        return self.atoms
    def push_atoms(self, atoms):
        for a in atoms:
            self.push_atom(a)
        return self.atoms
    def pop_atom(self, atom_index):
        txt = '\n'.join([self.atoms[x].info_txt() for x in sorted(atom_index,reverse=True)])
        [self.atoms.pop(x) for x in sorted(atom_index, reverse=True)]
        return txt
    def renumber_atom_id(self):
        for i,atom in enumerate(self.atoms):
            atom.atom_id = i+1
        return self.atoms
    def reverse_res_num(self):
        renumber_res_chain_id(False)
        resid_b = self.atoms[0].res_num
        resid_e = self.atoms[-1].res_num
        for atom in self.atoms:
            atom.res_num = resid_e - atom.res_num - resid_b
        return
    def renumber_res_chain_id(self, ch = True, res=True):
        if ch: self.renumber_chain_id()
        if res: self.renumber_res_num()
        return 0
    def renumber_res_num(self):
        self.set_residues_from_atom_info()
        for res_index, res in self.residues.items():
            for atom_id in res.atom_indice:
                self.atoms[atom_id].res_num = res_index
        return 0
    def renumber_chain_id(self):
        for chain_id, chain in self.chains.items():
            for atom in chain.atom_indice:
                atom.chain_id = chain_id
        return 0
    def push_atoms_renumber(self, atoms):
        for atom in atoms:
            self.push_atom(atom)
        self.renumber_res_chain_id(ch=False)
    def pop_residue_with_atom_id(self, atom_id):
        atom_index = self.get_atom_index_by_atom_id(atom_id)
        res_nums = set([self.atoms[idx].res_num for idx in atom_index])
        atom_index_pop = self.get_atom_index_by_res_num(res_nums)
        txt = self.pop_atom(atom_index_pop)
        txt += "\n" + str(len(atom_index_pop)) + " atoms were deleted\n"
        ##self.renumber_atom_id()
        return txt

    ### get atom_index ####
    def get_atom_index_by_atom_id(self, ids):
        d = {}
        for i,atom in enumerate(self.atoms):
            if atom.atom_id in ids:
                d[atom.atom_id] = i
        ret = []
        for aid in ids:  ret.append(d[aid])
        return ret
    def get_atom_index_by_atom_name(self, name, subset=[]):
        if len(subset) == 0:
            subset = range(0, len(self.atoms))
        ret = set()
        for index in subset:
            ##print "[" + self.atoms[index].atom_name + "] vs [" + name + "]"
            if self.atoms[index].atom_name == name:
                ret.add(index)
        return list(ret)
    #def get_atom_index_by_atom_id(self, atom_ids):
    #    atom_index = []
    #    for i,atom in enumerate(self.atoms):
    #        if atom.atom_id in atom_ids:
    #            atom_index.append(i)
    #    return atom_index
    def get_atom_index_by_res_num(self, res_nums):
        atom_index = []
        for i,atom in enumerate(self.atoms):
            if atom.res_num in res_nums:
                atom_index.append(i)
        return atom_index
    def get_atom_index_by_res_name(self, res_names):
        atom_index = []
        for i,atom in enumerate(self.atoms):
            if atom.res_name in res_names:
                atom_index.append(i)
        return atom_index

    def get_atom_index_by_res_num_and_atom_name(self, queries):
        atom_index = [-1] * len(queries)
        for i,atom in enumerate(self.atoms):
            for j,(q_res_num,q_atom_name) in enumerate(queries):
                if atom.res_num == q_res_num and \
                        atom.atom_name == q_atom_name:
                    atom_index[j] = i
        return atom_index
    def get_res_num_names(self):
        res_num_name = {}
        for atom in self.atoms:
            res_num_name[atom.res_num] = atom.res_name
        return res_num_name

    #### setter 

    def set_model_id(self, model_id):
        self.model_id = model_id
        return self.model_id
    def set_chain_id(self, chain_id):
        for atom in self.atoms:
            atom.chain_id = chain_id
        return
    def change_res_num_name(self, bef_id, aft_id, aft_name):
        for atom in self.atoms:
            if atom.res_num == bef_id:
                atom.res_num = aft_id
                atom.res_name = aft_name
        return 
    #### calc geometric features

    def get_center_self(self, index=[]):
        ##return self.get_center(self.atoms)
        ##def get_center(atoms):
        if len(index) == 0: index = range(0,len(self.atoms))
        crd = numpy.array([0.0, 0.0, 0.0])
        for i in index: crd += self.atoms[i].crd
        return crd/float(len(self.atoms))
    def get_circumscribed_box(self):
        box_min = numpy.array([1.0e10,   1.0e10,  1.0e10])
        box_max = numpy.array([-1.0e10, -1.0e10, -1.0e10])
        for atom in self.atoms:        
            if atom.crd[0] < box_min[0]: box_min[0] = atom.crd[0]
            if atom.crd[1] < box_min[1]: box_min[1] = atom.crd[1]
            if atom.crd[2] < box_min[2]: box_min[2] = atom.crd[2]
            if atom.crd[0] > box_max[0]: box_max[0] = atom.crd[0]
            if atom.crd[1] > box_max[1]: box_max[1] = atom.crd[1]
            if atom.crd[2] > box_max[2]: box_max[2] = atom.crd[2]
        return box_min,box_max

    #### modification geometries

    def translate(self,trans):
        for atom in self.atoms:
            atom.crd += trans
        #[atom.crd += trans for atom in self.atoms]
        return

    def rotate(self,rot):
        for atom in self.atoms:
            atom.crd = numpy.dot(rot,atom.crd)
        return
    def get_crd_matrix(self, index=[]):
        if len(index) == 0: index = range(0,len(self.atoms))
        return numpy.copy(numpy.matrix([self.atoms[i].crd for i in index]))
    def set_crd_matrix(self, mtx, index=[]):
        if len(index)==0: index=range(0,len(mtx))
        for j,i in enumerate(index):
            self.atoms[i].crd = mtx[j]
        return 

    def get_fsa(self):
        chains = {}
        for atom in self.atoms:
            if not atom.chain_id in chains:
                chains[atom.chain_id] = ""
            if atom.atom_name == "CA":
                chains[atom.chain_id] += kk.AA_3_1[atom.res_name]
            elif atom.atom_name == "C5'":
                chains[atom.chain_id] += kk.NA_3_1[atom.res_name]
                
        fasta = ""
        for chid, seq in chains.items():
            text =  ">CHAIN_" + str(chid) + "\n"
            text += "".join(seq) + "\n"
            fasta += text
        return fasta
     
    def checking_disconnect_of_residues(self, max_dist = 5.0, reassign_chain_id=False):
        assert len(self.chains.keys()) != 0
        assert len(self.residues.keys()) != 0

        ## disconnect[chain_id1][res_index1] = (flg_dist, flg_id)
        ## flg_dist (0 or 1): checking O-N distance 
        ## flg_id (0 or 1): checking res_num or res_name
        disconnect = {}

        def errout(atom_name, res):
            errmsg = "Error: Find unreasonable number of "+ atom_name +" atoms in residue "
            errmsg += str(res.res_num) + ":\n"
            for atom_index in res.atom_indice:
                errmsg += self.atoms[atom_index].info_txt() + "\n"
            sys.stderr.write(errmsg)
            sys.exit(1)

        for chain_id in sorted(self.chains.keys()):
            #print "CHAIN " + chain_id
            disconnect[chain_id] = {}
            a_name1 = "O"
            a_name2 = "N"
                
            if self.chains[chain_id].chain_type == kk.CH_DNA:
                a_name1 = "O3'"
                a_name2 = "P"
            elif not self.chains[chain_id].chain_type == kk.CH_PEPTIDE:
                if DEBUG: print "chain " + chain_id + " has undefined type"
                continue
            
            for res_index1 in self.chains[chain_id].res_indice:
                res_index2 = res_index1 + 1
                if not res_index1 in self.chains[chain_id].res_indice: continue
                if not res_index2 in self.chains[chain_id].res_indice:
                    #print "res " + str(res_index2) + " is not in chain " + chain_id + ":"
                    #print self.chains[chain_id].res_indice
                    disconnect[chain_id][res_index1] = (False, False)
                    continue
                if self.residues[res_index1].res_name == "NME":
                    disconnect[chain_id][res_index1] = (False, True)
                    continue
                if self.residues[res_index1].res_name == "ACE":
                    a_name1 = "C"
                res1a = self.get_atom_index_by_atom_name(a_name1, self.residues[res_index1].atom_indice)
                res2a = self.get_atom_index_by_atom_name(a_name2, self.residues[res_index2].atom_indice)
                try: assert len(res1a) == 1
                except:
                    errout(a_name1, self.residues[res_index1])
                if len(res2a) == 0:
                    disconnect[chain_id][res_index1] = (False, True)
                    continue
                try: assert len(res2a) == 1
                except:
                    errout(a_name2, self.residues[res_index2])

                res1a = res1a[0]
                res2a = res2a[0]
                flg_dist = self.atoms[res1a].distTo(self.atoms[res2a]) >= max_dist
                ##DEBUG
                if flg_dist:
                    print "distance : "+ str(res1a) +"-"+str(res2a) + " " + str(res_index1) + "-" + str(res_index2) + "  " + str(self.atoms[res1a].distTo(self.atoms[res2a]))
                flg_num = self.atoms[res1a].res_num != self.atoms[res2a].res_num - 1
                if flg_dist or flg_num:
                    disconnect[chain_id][res_index1] = (flg_dist, flg_num)
        if DEBUG:
            for chain_id, res_info in disconnect.items():
                for res_index, flgs in res_info.items():
                    if flgs[0] or flgs[1]:
                        line = str(chain_id) + ":" + str(self.residues[res_index].res_num) + "("+str(res_index) + ") is disconnected: "
                        if flgs[0]: line += " distance"
                        if flgs[1]: line += " residue_id"
                        #print line

        if reassign_chain_id:
            self.reassign_chain_id_to_atoms(disconnect)
        return 0

    def reassign_chain_id_to_atoms(self,disconnect):
        offset = 0
        res_index_new_chain = {}
        term_res = []
        #print disconnect
        for chain_id in sorted(disconnect.keys()):
            res_indice = disconnect[chain_id]
            for res_index, flgs in res_indice.items():
                term_res.append(res_index)
        term_res = sorted(term_res)

        term_res_chid = {}
        pointer_chain_id = 'A'
        for res in term_res:
            term_res_chid[res] = pointer_chain_id
            pointer_chain_id = chr(ord(pointer_chain_id)+1)

        #print term_res_chid
        for res_index, residue in self.residues.items():
            new_ch_id = "X"
            for i, cur_res in enumerate(term_res):
                if res_index <= cur_res:
                    new_ch_id = term_res_chid[cur_res]
                    break
            for atom_index in residue.atom_indice:
                self.atoms[atom_index].chain_id = new_ch_id

        #for chain_id in sorted(disconnect.keys()):
        #    res_indice = disconnect[chain_id]
        #    for res_index, flgs in res_indice.items():
        #        pointer_chain_id = chr(ord(chain_id)+1)
        #        if flgs[0] or flgs[1]: offset += 1
        #        res_index_new_chain[res_index] = pointer_chain_id
        #for res_index, residue in self.residues.items():
        #    for atom_index in residue.atom_indice:
        #        self.atoms[atom_index].chain_id = res_index_new_chain[res_index]
        return 0
    
    def swap_res_num_auth(self):
        for res_num, res in self.residues.items():
            tmp = res.res_num
            res.res_num = res.res_num_auth
            res.res_num_auth = tmp
        for atom in self.atoms:
            tmp = atom.res_num
            atom.res_num = atom.res_num_auth
            atom.res_num_auth = tmp
        return 
    def set_res_num_auth_from_model(self, auth_model, skip_ids = None):
        resid_conv = {}
        res_num_auth = -1
        for res_num, res in self.residues.items():
            if skip_ids and res_num in skip_ids: continue
            res_num_auth += 1

            if not res_num_auth in auth_model.residues: break
                
            if (res.res_name[0:3] == auth_model.residues[res_num_auth].res_name[0:3]) or \
                    (kk.RESNAMES[res.res_name] == kk.RESNAMES[auth_model.residues[res_num_auth].res_name]):
                self.residues[res_num].res_num_auth = auth_model.residues[res_num_auth].res_num
                self.residues[res_num].res_name = auth_model.residues[res_num_auth].res_name
                #print str(res_num) + ":" + res.res_name + " - " + str(res_num_auth) + ":" + auth_model.residues[res_num_auth].res_name + ":" + str(auth_model.residues[res_num_auth].res_num)
                self.residues[res_num].chain_id = auth_model.residues[res_num_auth].chain_id
                
        for i, atom in enumerate(self.atoms):
            self.atoms[i].res_num_auth = self.residues[atom.res_index].res_num_auth
            self.atoms[i].chain_id = self.residues[atom.res_index].chain_id

        self.reset_chains_from_atom_info()
        return
    def set_water_chain(self, wat_chain_id="X"):
        print "WATER"
        for i, atom in enumerate(self.atoms):
            if kk.RESNAMES[atom.res_name] == "HOH":
                self.chains[atom.chain_id].res_indice.remove(atom.res_index)
                self.chains[atom.chain_id].atom_indice.remove(i)
                if not wat_chain_id in self.chains:
                    self.chains[wat_chain_id] = Chain(wat_chain_id)
                self.chains[wat_chain_id].push_res_index(atom.res_index)
                self.chains[wat_chain_id].push_atom_index(i)
                atom.chain_id = wat_chain_id
                print atom.info_txt()
        return 
    def insert_atoms_res(res_index, atoms):

        return
