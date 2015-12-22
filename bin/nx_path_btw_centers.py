#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import networkx as nx
import re
import numpy as np
from optparse import OptionParser
import sys
import os
import sqlite3

sys.path.append(os.environ["HOME"]+"/local/kktools/kkio")
import kktable
sys.path.append(os.environ["HOME"]+"/local/kktools/network")
from nx_centrality import read_headers

def define_options():
    p = OptionParser()
    p.add_option('-i', dest='fn_table',
                 help="kktable file describing elements")
    p.add_option('--i-rel', dest='fn_table_rel',
                 help="kktable file describing relations")
    p.add_option('-o', dest='fn_out',
                 help="Output file")
    p.add_option('--o-edgelist', dest='fn_o_edge_list',
                 help="Output file")
    p.add_option('--col', dest='col_centrality',
                 help="Column name for centrality")
    p.add_option('--min-cent', dest='min_centrality',
                 type="float",
                 help="Minimum centrality")

    p.add_option('--top-atoms', dest='top_atoms',
                 type="float",
                 help="Parcentage of number of atoms treated as center atoms")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = define_options() 

    headers, initial_ch = read_headers(opts.fn_table_rel)
    print headers
    print "read_edgelist"
    graph = nx.read_edgelist(opts.fn_table_rel,
                             comments=initial_ch,
                             nodetype=headers[0][1],
                             data=headers[2:])

    elems, type_def = kktable.TableReader(opts.fn_table).read_binrel_elements("atom_id")
    print "table: " + str(len(elems.keys()))
    min_cent = 0
    if opts.min_centrality:
        min_cent = opts.min_centrality
    if opts.top_atoms:
        cent_vals = []
        for atom_id, keyval in elems.items():
            cent_vals.append(keyval[opts.col_centrality])
            cent_vals.sort(reverse=True)
        n_atoms = int(opts.top_atoms * len(cent_vals))
        min_cent = cent_vals[n_atoms-1]

    print "Minimum centrality: " + str(min_cent)
    center_atoms = []
    for atom_id, keyval in elems.items():
        if keyval[opts.col_centrality] > min_cent:
            center_atoms.append(atom_id)

    print "Center atoms: " + str(len(center_atoms)) + " / " + str(len(elems.keys()))

    f_out = None
    if opts.fn_out:
        f_out = open(opts.fn_out, "w")
    
    edge_set = set()
    for atom_id1 in center_atoms:
        for atom_id2 in center_atoms:
            if atom_id1 == atom_id2: break
            path = None
            try:
                path = nx.shortest_path(graph, atom_id1, atom_id2)
            except nx.NetworkXNoPath:
                continue
            #for path in paths:
            for i in range(len(path)-1):
                edge = (path[i], path[i+1])
                if path[i] > path[i+1]:
                    edge = (path[i+1], path[i])                    
                edge_set.add(edge)
            if f_out:
                f_out.write(" ".join([str(x) for x in path]) + "\n")
    if f_out:  f_out.close()

    if opts.fn_o_edge_list:
        f_o_edge = open(opts.fn_o_edge_list,"w")
        for node1, node2 in edge_set:
            f_o_edge.write(str(node1)+"\t"+str(node2)+"\t1\n")
        f_o_edge.close()
    return 

if __name__ == "__main__":
    _main()
