#!/usr/bin/python2.7

import kkkit
import kkbinrel
import re
import numpy as np
import sqlite3

class TableReadManager(object):
    def __init__(self, fn_elem, fn_rel, elem_key, rel_key_source=None, rel_key_target=None):
        self.fn_elem = fn_elem
        self.fn_rel = fn_rel
        self.elem_key = elem_key
        self.key_source = rel_key_source
        self.key_target = rel_key_target
        return 
    def read_table_matrix(self):
        binrel = kkbinrel.BinaryRelation(elem_key=self.elem_key)
        elem, type_def = TableReader(self.fn_elem).read_binrel_elements(self.elem_key)
        rel, type_def_rel = TableReader(self.fn_rel).read_binrel_relation_matrix(self.elem_key)
        binrel.elements = elem
        binrel.relations = rel
        binrel.push_type_defs(type_def)
        binrel.push_type_defs(type_def_rel)
        return binrel
    def read_tables(self, sym):
        binrel = kkbinrel.BinaryRelation(symmetry=sym, elem_key=self.elem_key)
        elem, type_def = TableReader(self.fn_elem).read_binrel_elements(self.elem_key)
        rel, type_def_rel = TableReader(self.fn_rel).read_binrel_relations(key_source=self.key_source, key_target=self.key_target)
        binrel.elements = elem
        binrel.relations = rel
        binrel.push_type_defs(type_def)
        binrel.push_type_defs(type_def_rel)
        return binrel

class TableReaderDB():
    def __init__(self, fn):

        return
    def read_binrel_elements(self, id_key = ""):
        self.con = sqlite3.connect(fn)
        sql = ""
        self.con.close()
        return 

class TableReader(kkkit.FileI):
    def __init__(self, fn):
        super(TableReader, self).__init__(fn)
        self.na_string = "N/A"
    def read_binrel_elements(self, id_key = ""):
        self.open()
        line = self.f.readline()
        headers, type_def = self.read_header(line)

        keys, type_def = self.read_header(line)
        binrel = kkbinrel.BinaryRelation()
        binrel.push_type_defs(type_def)
        
        if id_key == "":
            id_key = keys[0]
        id_col_num = headers.index(id_key)
        if id_col_num == -1: id_col_num = 0
        for line in self.f:
            terms = re.compile("\s+").split(line.strip())
            elem = {}
            for i,val in enumerate(terms):
                elem[keys[i]] = val
            binrel.push_element(elem[headers[id_col_num]], elem)
        self.close()
        return binrel.elements, binrel.type_def

    def read_binrel_relation_matrix(self, key):
        self.open()
        line = self.f.readline()
        headers = re.compile("\s+").split(line.strip())
        binrel = kkbinrel.BinaryRelation()
        for line in self.f:
            terms = re.compile("\s+").split(line.strip())
            rel = {}
            for i,val in enumerate(terms[1:]):
                if val == self.na_string: continue
                keyval = {key:val}
                binrel.push_relation(terms[0], headers[i], keyval)
        self.close()
        return binrel.relations

    def read_header(self, line):
        keys = []
        type_def = {}
        headers = re.compile("\s+").split(line.strip())
        for header in headers:
            terms = re.compile("\.").split(header)
            if len(terms) > 1 and \
                    terms[-1] in kkbinrel.BinaryRelation.types:
                key = ".".join(terms[:-1])
                type_def[key] = terms[-1]
                keys.append(key)
            else: keys.append(header)
        return keys, type_def
    def read_binrel_relations(self, 
                              abs_min_key="", abs_min_val=0.0,
                              max_key="", max_val=0.0,
                              minmax_filter_and=False,
                              key_source=None, key_target=None,
                              larger_priority_key=None, prio_func=np.fabs):
        ## abs_min_key, abs_min_val are filter setting
        ## Only when the absolute value of attribute named abs_min_key
        ## is larger than abs_min_val, the relationship is read
        ## These options are used to filt out relations with
        ## low correlation coefficient.

        ## c_source, c_target
        ## specifies index of column meaning 
        ## the source and target nodes

        ## larger_priority_key
        ## should be specified as key name
        ## If there are some edges with identical source and target,
        ## edge with the largest values in this key
        ## is used


        self.open()
        line = self.f.readline()
        keys, type_def = self.read_header(line)
        binrel = kkbinrel.BinaryRelation()
        binrel.push_type_defs(type_def)

        n_rel = 0
        c_source = -1
        c_target = -1
        if key_source and key_target:
            flg = False
            for i, key in enumerate(keys):
                if key == key_source: c_source = i
                elif key == key_target: c_target = i
            if c_source==-1:
                print "column ["+key_source+"] is not defined" 
                flg = True
            if c_target==-1:
                print "column ["+key_target+"] is not defined"
                flg = True
            if flg:
                return
        else:
            c_source = 0
            c_target = 1
            
        for line in self.f:
            terms = re.compile("\s+").split(line.strip())
            keyval = {}
            for i,val in enumerate(terms):
                if i != c_source and i != c_target:
                    keyval[keys[i]] = val
            flg_read = True

            flg_min = True
            if abs_min_key in keyval and \
                    np.fabs(float(keyval[abs_min_key])) < abs_min_val:
                flg_min = False
            flg_max = True
            if max_key in keyval and \
                    float(keyval[max_key]) > max_val:
                flg_max = False
            if minmax_filter_and:
                flg_read = flg_min & flg_max
            else:
                flg_read = flg_min | flg_max

            pair = (terms[c_source], terms[c_target])
            if larger_priority_key and larger_priority_key in keyval and \
                    pair in binrel.relations and \
                    prio_func(binrel.relations[pair][larger_priority_key]) >= \
                    prio_func(keyval[larger_priority_key]):
                flg_read = False

            if flg_read:
                binrel.push_relation(terms[c_source], terms[c_target], keyval)
                n_rel += 1
        print str(n_rel) + " relations were read."
        self.close()
        return binrel.relations, binrel.type_def


class TableWriter(kkkit.FileO):
    def __init__(self, fn):
        super(TableWriter, self).__init__(fn)
        self.na_string = "N/A"
    def write_binrel_elements(self, binrel):
        keys = set()
        self.open()
        for elem_id, elem_keyval in binrel.elements.items():
            for key, val in elem_keyval.items():
                keys.add(key)
        keys = list(keys)
        headers = [binrel.element_key + "." + binrel.element_key_type]
        for key in keys:
            if key != binrel.element_key:
                headers.append(key + "." + binrel.type_def[key])
            
        self.f.write("\t".join(headers)+"\n")
        for elem_id, elem_keyval in binrel.elements.items():
            terms = [str(elem_id)]
            for key in keys:
                if key == binrel.element_key: continue
                v = self.na_string
                if key in elem_keyval:
                    v = elem_keyval[key]
                terms.append(str(v))
            self.f.write("\t".join(terms) + "\n")
        self.close()
        return 

    def write_binrel_relations(self, binrel):
        keys = set()
        self.open()
        for rel_id, rel_keyval in binrel.relations.items():
            for key, val in rel_keyval.items():
                keys.add(key)
        keys = list(keys)
        headers = [binrel.relation_key1 + "." + binrel.element_key_type,
                   binrel.relation_key2 + "." + binrel.element_key_type]
        for key in keys:
            headers.append(key + "." + binrel.type_def[key])
        
        self.f.write("\t".join(headers)+"\n")
        for rel_id, rel_keyval in binrel.relations.items():
            terms = [str(rel_id[0]), str(rel_id[1])]
            for key in keys:
                v = self.na_string
                if key in rel_keyval: v = rel_keyval[key]
                terms.append(str(v))
            self.f.write("\t".join(terms) + "\n")
        self.close()
        return 
        
    def write_binrel_relation_matrix(self, binrel, key):
        self.open()
        elems = list(binrel.elements.keys())
        line = "\t".join([str(x) for x in elems])
        self.f.write(line+"\n")
        for elem1 in elems:
            line = str(elem1)
            for elem2 in elems:
                elem_pair = (elem1, elem2)
                if elem1 > elem2:
                    elem_pair = (elem2, elem1)
                if (elem_pair in binrel.relations) and \
                        key in binrel.relations[elem_pair]:
                    line += "\t" + str(binrel.relations[elem_pair][key])
                else:
                    line += "\t" + self.na_string
            line += "\n"
            self.f.write(line)
        self.close()
        return
    
