#!/usr/bin/python2.7

import kkkit
import kkstruct
import re

class PDBWriter(kkkit.FileO):
    def __init__(self, fn):
        super(PDBWriter, self).__init__(fn)
    def write_model(self, model, flg_presto=False):
        self.open()
        if model.title != "":
            txt = "TITLE " + model.title
            self.f.write(txt+'\n')
        if model.pbc_box != [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]:
            txt = "CRYST1%9.2f%9.2f%9.2f%7.2f%7.2f%7.2f"%(
                model.pbc_box[0],
                model.pbc_box[1],
                model.pbc_box[2],
                model.pbc_box[3],
                model.pbc_box[4],
                model.pbc_box[5]
                )
            txt += '  P 1           1           '
            self.f.write(txt+'\n')            
        for atom in model.atoms:
            atom_name = atom.atom_name
            if len(atom_name) < 4: atom_name = " " + atom_name
            assert(len(atom.chain_id)==1)
            txt = "ATOM"
            if flg_presto:
                if atom.atom_id <= 99999:
                    txt += "  "
                elif atom.atom_id <= 999999:
                    txt += " "
                txt += "%5d"%atom.atom_id #6:11
            else:
                txt += "  "
                aid = atom.atom_id
                if atom.atom_id > 99999:
                    aid = atom.atom_id%100000
                txt += "%5d"%aid #6:11
            if atom.res_num > 9999:
                atom.res_num -= 10000
                atom.chain_id = chr(ord(atom.chain_id[0])+1)

            txt += " "
            txt += "%-4s"%atom_name #12:16
            txt += atom.alt_loc[0]
            txt += "%-4s"%atom.res_name #17:20
            txt += atom.chain_id[0] #21
            txt += "%4d    "%atom.res_num #22:26
            txt += "%8.3f"%atom.crd[0] #30:38
            txt += "%8.3f"%atom.crd[1] #38:46
            txt += "%8.3f"%atom.crd[2] #46:54
            txt += "%6.2f"%atom.ocp #54:60
            txt += "%6.2f      "%atom.tf #60:66
            txt += "%4s"%atom.seg_id #72:76
            txt += "%2s"%atom.elem #76:78
            self.f.write(txt+'\n')
        self.close()
        return

class PDBReader(kkkit.FileI):
    def __init__(self, fn):
        super(PDBReader, self).__init__(fn)
    def read_model(self, flg_ignore_exception=False):
        self.open()
        model = kkstruct.Model()
        model_id = 0
        atom_id = -1
        line = self.f.readline()
        while line:
            atom_id_pdb = 0
            atom_name = ""
            res_name = ""
            chain_id = ""
            res_num = 0
            x = 0.0
            y = 0.0
            z = 0.0
            ocp = 0.0
            tf = 0.0
            seq_id = ""
            elem = ""
            chg = 0.0
            header = line[0:6]
            if header[0:5] == 'MODEL':
                model_id = int(line[6:])
                atom_id = -1
                model.set_model_id(model_id)
            elif header[0:5] == "TITLE":
                model.title = line[6:].strip()
            elif header[0:6] == "CRYST1":
                model.pbc_box = [float(line[6:15]),
                                 float(line[15:24]),
                                 float(line[24:33]),
                                 float(line[33:40]),
                                 float(line[40:47]),
                                 float(line[47:54])]
            elif header[0:4] == 'ATOM' or header == 'HETATM':
                atom_id += 1
                ma1 = re.compile("\d+").match(header[4])
                ma2 = re.compile("\d+").match(header[5])
                beg_atom_id = 6
                if ma1 and ma2: beg_atom_id = 4
                if not ma1 and ma2: beg_atom_id = 5
                try: atom_id_pdb = int(line[beg_atom_id:11])
                except: pass
                try: atom_name = line[12:16].strip()
                except: pass
                try: res_name = line[17:21].strip()
                except: pass
                try: chain_id = line[21]
                except: pass
                try: res_num = int(line[22:27])
                except: pass
                try: x = float(line[30:38])
                except: pass
                try: y = float(line[38:46])
                except: pass
                try: z = float(line[46:54])
                except: pass
                try: ocp = float(line[54:60])
                except: pass
                try: tf = float(line[60:66])
                except: pass
                try: seg_id = line[72:76]
                except: pass
                try: elem = line[76:78].strip()
                except: pass
                chg = 0.0
                #if len(line)>=80 and line[78:80] != "  ":
                #chg = float(line[78:80])
                model.push_atom(kkstruct.Atom(header,atom_id,atom_name,res_name,chain_id,res_num,x,y,z,ocp,tf,seg_id,elem,chg,atom_id_pdb))
            elif header == 'ENDMDL':
                break
            line = self.f.readline()
        self.close()
        model.set_residues_from_atom_info()
        model.set_chain_types()
        return model

