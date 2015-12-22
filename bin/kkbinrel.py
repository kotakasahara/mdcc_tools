#!/usr/bin/python2.7

import kkkit
import numpy
import kkdefine
import re
import sys

class BinaryRelation(object):
    types = {"int", "float", "string"}
    def __init__(self, symmetry=False, elem_key = "id"):
        """ Setting member variables
        """
        ## definition of tyep of each key
        self.type_def = {}

        ## symmetricity of matrix
        ## 0 or 1
        self.symmetry = symmetry

        ## metadata as a key-value dictionary
        self.meta = {}

        ## self.elements[id] = key-value dictionary
        self.elements = {}
        
        self.element_key = elem_key
        self.element_key_type = "int"

        ## self.relation[(elem_id1, elem_id2)] = key-value
        ## elem_id1 < elem_id2 for symmetry relationship
        self.relations = {}

        self.relation_key1 = self.element_key + "1"
        self.relation_key2 = self.element_key + "2"

        ## index of tables
        self.index_elem = {}
        self.index_rel = {}

        self.na_string = "N/A"
    def infer_type(self, val):
        if type(val) == type(1): return val, "int"
        elif type(val) != type("str"): return val, "float"
        tp = "string"
        m1 = re.compile("^-{0,1}\d+$").match(val)
        m2 = re.compile("^-{0,1}\d+\.\d+[Ee]{0,1}\d*$").match(val)
        if m1:
            val=int(val)
            tp = "int"
        elif m2:
            val = float(val)
            tp = "float"
        return val, tp
    def push_type_defs(self, defs):
        for key, define in defs.items():
            self.push_type_def(key,define)
        return 0
    def push_type_def(self, key, t):
        if t != "int" and t != "float" and t != "string":
            print "UNKNOWN TYPE : " + str(t)
            return 1
        self.type_def[key] = t
        if key == self.element_key_type:
            self.element_key_type = t            
        return 0

    def set_relations_from_matrix(self, matrix, elems, key, type_def,
                                  min_threshold=-1e100, max_threshold=1e100,
                                  flg_absolute = False):
        n_add = 0
        for i, row in enumerate(matrix):
            for j, elem in enumerate(row):
                val = elem
                if flg_absolute: val = numpy.fabs(val)
                if val >= min_threshold and val < max_threshold:
                    self.push_relation(elems[i], elems[j], {key:elem})
                    n_add += 1
        self.push_type_def(key, type_def)
        return n_add

    def make_index(self, table, key):
        index = {}
        for elem_id, entry in enumerate(table):
            index[entry[key]] = elem_id
        return index

    def set_type(self, val, tp):
        v = None
        try:
            if tp == "string":
                v = str(val)
            elif tp == "int":
                v = int(val)
            elif tp == "float":
                v = float(val)
        except:pass
        return v
    def set_types(self, keyval):
        for k, v in keyval.items():
            if v==self.na_string: continue
            if k in self.type_def:
                v = self.set_type(v, self.type_def[k])
                if v != None:
                    keyval[k] = v
            else:
                v, tp = self.infer_type(v)
                #self.push_type_def(k, tp)
                self.type_def[k] = tp
                keyval[k] = self.set_type(v, tp)
        return keyval
    def push_element(self, i, new_elem):
        i = self.set_type(i, self.element_key_type)        
        new_elem = self.set_types(new_elem)
        ## new_elem is key-value dictionary
        if i in self.elements:
            for k,v in new_elem.items():
                self.elements[i][k] = v
        else:
            self.elements[i] = new_elem
        return len(self.elements.keys())

    def push_relation(self, i, j, keyval):
        i = self.set_type(i, self.element_key_type)
        j = self.set_type(j, self.element_key_type)
        keyval = self.set_types(keyval)

        ## keyval is key-value dictionary
        if (i,j) in self.relations:
            for k,v in keyval.items():
                self.relations[(i,j)][k] = v
        else:
            self.relations[(i,j)] = keyval
        return 0

    def get_matrix(self, key, default=0.0, diagonal=1.0):
        mtx = []
        for i,elem_i in self.elements.items():
            row = []
            for j,elem_j in self.elements.items():
                pair = (i, j)
                if self.symmetry and j<i: pair = (j,i)

                val = default
                if i==j: val = diagonal

                if pair in self.relations:
                    val = self.relations[pair][key]
                row.append(val)
            mtx.append(numpy.array(row))
        return numpy.array(mtx)

    def get_fsa_from_elem_res(self):
        """
        When each element corresponds one residue,
        generate fasta file for each chain_id
        """
        chains = {}
        for i, elem in self.elements.items():
            if not elem["chain_id"] in chains:
                chains[elem["chain_id"]] = ""
            if elem["atom_name"] == "CA":
                chains[elem["chain_id"]] += kkdefine.AA_3_1[elem["res_name"]]
            elif elem["atom_name"] == "C5'":
                chains[elem["chain_id"]] += kkdefine.NA_3_1[elem["res_name"]]
        fasta = ""
        for chid, seq in chains.items():
            text =  ">CHAIN_" + str(chid) + "\n"
            text += "".join(seq) + "\n"
            fasta += text
        return fasta

    def union_elem_bykey(self, oth, ekey, prio_arg=False):
        self.union_dict(self.meta, oth.meta, prio_arg)
       ## Elements

        for oth_elem_id, oth_elem in oth.elements.items():
            for elem_id, elem in self.elements.items():
                if elem[ekey] == oth_elem[ekey]:
                    self.union_dict(elem, oth_elem, prio_arg)
        ###Type def
        self.union_dict(self.type_def, oth.type_def, prio_arg)
        return 

    def union(self, oth, prio_arg=False):
        """
        Merging data in two BinaryRelation class.
        prio_arg=True means values in self is 
        overwritten by oth's value
        when they are conflicts.
        Alignment of elements is performed by key_to_elem_id
        """
        self.union_dict(self.meta, oth.meta, prio_arg)

        if self.element_key != oth.element_key:
            sys.stderr.write("Element keys are in consistent:")
            sys.stderr.write(self.element_key + " : " + oth.element_key)
            return None

        ## Elements
        for oth_elem_id, oth_elem in oth.elements.items():
            self.push_element(oth_elem_id, oth_elem)

        #### Relations
        print "LEN REL " + str(len(self.relations.keys()))
        for oth_rel_id, oth_rel in oth.relations.items():
            self.push_relation(oth_rel_id[0], oth_rel_id[1], oth_rel)
        
        ###Type def
        self.union_dict(self.type_def, oth.type_def, prio_arg)

        return 

    def union_dict(self, left, right, prio_arg=False):
        """
        this function called by self.union().
        argumetn left should be self member
        """
        #print left 
        #print right
        for key, val in right.items():
            left[key] = val
            #print key + " - " + str(val)
        return left
    def delete_rel_with_max_threshold(self, key, max_value):
        relations = {}
        for pair, keyval in self.relations.items():
            if key in keyval and keyval[key] < max_value:
                relations[pair] = keyval
        self.relations = relations
        return 
    def delete_rel_with_min_threshold(self, key, min_value):
        relations = {}
        for pair, keyval in self.relations.items():
            if key in keyval and keyval[key] >= min_value:
                relations[pair] = keyval
        self.relations = relations
        return 
