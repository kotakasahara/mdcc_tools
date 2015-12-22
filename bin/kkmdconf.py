#!/usr/bin/python2.7

import kkkit
import sys
import re

class MDConf(object):
    def __init__(self):
        self.description = ""
        self.gro_init_structure = ""
        self.gro_traj_edr = []
        self.gro_traj_trr = []
        self.gro_traj_set = []
        self.atom_selections = []
        self.selection_names = []
        self.segments = {}
        self.straddling = {}
        self.n_digit = 4
        self.gro_index_file = ""
        self.gro_tpr_file = ""
        self.psy_traj_crd = []
        self.psy_traj_vel = []
        self.psy_traj_set = []
        self.psy_init_structure = ""
        self.kk_traj_trans = ""
        self.pseudo_atoms = []
        
    def add_gro_traj_set(self, pref, begin, last):
        self.gro_traj_set.append([pref, begin, last])
    def add_gro_traj_file(self,pref):
        fn_edr = pref + ".edr"
        self.gro_traj_edr.append(fn_edr)        
        fn_trr = pref + ".trr"
        self.gro_traj_trr.append(fn_trr)
    def set_gro_traj_files(self):
        for info in self.gro_traj_set:
            st = "%0"+str(self.n_digit)+"d"
            for i in range(info[1],info[2]+1):
                run_id = st%i
                fn = info[0] + run_id
                self.add_gro_traj_file(fn)
    def add_psy_traj_set(self, pref, begin, last, suff_crd, suff_vel=""):
        self.psy_traj_set.append([pref, begin, last, suff_crd, suff_vel])
    def set_psy_traj_files(self):
        for info in self.psy_traj_set:
            st = "%0"+str(self.n_digit)+"d"
            for i in range(info[1],info[2]+1):
                suff_crd = info[3]
                suff_vel = info[4]
                run_id = st%i
                fn = info[0] + run_id
                self.psy_traj_crd.append(fn + "." + suff_crd)
                if suff_vel != "":
                    self.psy_traj_vel.append(fn + "." + suff_vel)
    def get_segment_bynum(self, num, default=None):
        for segname, seg in self.segments.items():
            if num >= seg[0] and num <= seg[1]:
                return segname, seg
        return default, (default,default,default)
    def set_straddling(self, seg_name, xyz):
        self.straddling[seg_name] = xyz
    def set_gro_index_file(self, fn):
        self.gro_index_file = fn
    def set_gro_tpr_file(self, fn):
        self.gro_tpr_file = fn
    def write_text(self):
        text = ""
        text += "--description--\n"
        text += self.description + "\n"
        text += "--description--\n"
        for fn in self.gro_traj_edr:
            text += fn + "\n"
        return text

class MDConfReader(kkkit.FileI):
    def __init__(self,fn):
        super(MDConfReader, self).__init__(fn)
    def read_conf(self):
        mdconf = MDConf()
        self.open()
        current_field = ""
        current_text = ""
        for orig_line in self.f:
            line = kkkit.eliminate_comment(orig_line).strip()
            if line.lower() == current_field:
                if current_field == "--description--":
                    mdconf.description = current_text
                current_field = ""
                current_text = ""
            elif line[0:2] == "--" and line[-2:] == "--":
                current_field = line.lower()
                
            elif current_field != "":
                current_text += line + "\n"
            else:
                terms = re.compile("\s+").split(line.strip())
                if len(line) == 0:
                    continue
                elif terms[0][0:2] != "--":
                    sys.stderr.write("Syntax error: Field tag must begin with '--':")
                    sys.stderr.write(line + "\n")
                elif terms[0] == "--gro-traj-set":
                    mdconf.add_gro_traj_set(terms[1],int(terms[2]),int(terms[3]))
                elif terms[0] == "--gro-traj":
                    mdconf.add_gro_traj_file(terms[1])
                elif terms[0] == "--gro-index":
                    mdconf.set_gro_index_file(terms[1])
                elif terms[0] == "--gro-tpr":
                    mdconf.set_gro_tpr_file(terms[1])
                elif terms[0] == "--gro-initial-structure":
                    mdconf.gro_init_structure = terms[1]
                elif terms[0] == "--atom-selection":
                    mdconf.atom_selections.append(" ".join(terms[1:]))
                elif terms[0] == "--selection-name":
                    mdconf.selection_names.append(" ".join(terms[1:]))
                elif terms[0] == "--segment":
                    seg = (int(terms[2]), int(terms[3]), terms[4])
                    mdconf.segments[terms[1]] = seg
                elif terms[0] == "--straddling":
                    mdconf.set_straddling(terms[1], [int(x) for x in terms[2:5]])
                elif terms[0] == "--n-digit-run-id":
                    mdconf.n_digit = int(terms[1])
                elif terms[0] == "--psy-traj-set":
                    if len(terms)==5:
                        mdconf.add_psy_traj_set(terms[1],int(terms[2]),int(terms[3]),terms[4])
                    elif len(terms) >= 6:
                        mdconf.add_psy_traj_set(terms[1],int(terms[2]),int(terms[3]),terms[4],terms[5])
                elif terms[0] == "--psy-initial-structure":
                    mdconf.psy_init_structure = terms[1]
                elif terms[0] == "--kk-traj-trans":
                    mdconf.kk_traj_trans = terms[1]
                elif terms[0] == "--pseudo-atom":
                    mdconf.pseudo_atoms.append([float(x) for x in terms[1:4]])
                else:
                    sys.stderr.write("Unknown tag:")
                    sys.stderr.write(terms[0])
        self.close()
        mdconf.set_gro_traj_files()
        return mdconf
