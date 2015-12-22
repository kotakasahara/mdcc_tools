#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import networkx as nx
import re
import numpy as np
from optparse import OptionParser
import sys
import os
import kktable 
import kkbinrel

def define_options():
    p = OptionParser()
    p.add_option('-i', dest='fn_table',
                 help="kktable file describing relations")
    p.add_option('--i-elem', dest='fn_elements',
                 help="kktable file describing elements")
    p.add_option('--key-elem', dest='key_element',
                 default="node_id",
                 help="key of elements")
    p.add_option('--key-src', dest='key_src',
                 default="node_id1",
                 help="key of elements")
    p.add_option('--key-dst', dest='key_dst',
                 default="node_id2",
                 help="key of elements")
    p.add_option('-o', dest='fn_out',
                 help="Output file")

    p.add_option('--o-group', dest='fn_out_group',
                 help="Output file for grouped nodes")
    p.add_option('--key-group', dest='key_group',
                 help="key of node gouprs")
    p.add_option('--i-group', dest='fn_groups',
                 help="gropu tables")

    p.add_option('--o-header', dest='flg_out_header',
                 action="store_true",
                 help="A header line will be added to the output file.")

    p.add_option('--dgr', dest='flg_degree',
                 action="store_true",
                 help="Calculate degree centrality")
    p.add_option('--cls', dest='flg_closeness',
                 action="store_true",
                 help="Calculate closeness centrality")
    p.add_option('--btw', dest='flg_betweenness',
                 action="store_true",
                 help="Calculate betweenness centrality")
    p.add_option('--inf', dest='flg_information',
                 action="store_true",
                 help="Calculate information centrality")
    p.add_option('--rwb', dest='flg_randomwalk',
                 action="store_true",
                 help="Calculate randomwalk betweenness centrality")
    p.add_option('--egv', dest='flg_eigenvector',
                 action="store_true",
                 help="Calculate eigenvector centrality")
    p.add_option('--cmn', dest='flg_communicability',
                 action="store_true",
                 help="Calculate communicability centrality")
    p.add_option('--lod', dest='flg_load',
                 action="store_true",
                 help="Calculate load centrality")

    p.add_option('--key-dist', dest='key_distance',
                 help="key to distances of edges")
    p.add_option('--key-weight', dest='key_weight',
                 help="key to weight of edges")

    p.add_option('--abs-key', dest='abs_key',
                 action="append",
                 help="key to be converted into their absolute values")
    p.add_option('--rcp-key', dest='rcp_key',
                 action="append",
                 help="key to be converted into their reciprocal values")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = define_options() 
    headers, initial_ch = read_headers(opts.fn_table)
    print headers
    print "read_edgelist"
    
    elems = {}
    elem_types = {}
    if opts.fn_elements:
       elems, elem_types = kktable.TableReader(opts.fn_elements).read_binrel_elements(opts.key_element)
    print "Elements:" + str(len(elems.keys()))

    graph = nx.read_edgelist(opts.fn_table,
                             comments=initial_ch,
                             nodetype=headers[0][1],
                             data=headers[2:])
    if opts.abs_key:
        print "add_field abs"
        graph = add_field(graph, opts.abs_key, np.fabs, "abs_")
    if opts.rcp_key:
        print "add_field rcp"
        graph = add_field(graph, opts.rcp_key, rcp, "rcp_")
    
    centralities = {}
    cent_keys = []

    if opts.flg_degree:
        print "degree centrality"
        deg = nx.degree_centrality(graph)
        centralities["degree"] = deg
        cent_keys.append("degree")
    if opts.flg_closeness:
        print "closeness centrality"
        if opts.key_distance: print "distance: " + opts.key_distance
        cls = nx.closeness_centrality(graph, distance=opts.key_distance)
        centralities["closeness"] = cls
        cent_keys.append("closeness")
    if opts.flg_betweenness:
        print "betweenness centrality"
        if opts.key_weight: print "weight: " + opts.key_weight
        btw = nx.betweenness_centrality(graph, weight=opts.key_weight)
        centralities["betweenness"] = btw
        cent_keys.append("betweenness")
    if opts.flg_information:
        print "information centrality (current flow closeness )"
        if opts.key_weight:
            print "weight: " + opts.key_weight
            try:
                inf = nx.current_flow_closeness_centrality(graph, weight=opts.key_weight)
                centralities["information"] = inf
                cent_keys.append("information")
            except:
                print "ERROR in information centrality"
        else:
            print "--key-weight is required for calculating information centrality"
    if opts.flg_randomwalk:
        print "random walk betweenness centrality ( current flow_betweenness )"
        if opts.key_weight:
            print "weight: " + opts.key_weight
            try:
                rwb = nx.current_flow_betweenness_centrality(graph, weight=opts.key_weight)
                centralities["randomwalk_betweenness"] = rwb
                cent_keys.append("randomwalk_betweenness")
            except:
                print "ERROR in random walk betweenness centrality"
        else:
            print "--key-weight is required for calculating random-walk betweenness centrality"
    if opts.flg_eigenvector:
        print "eigenvector centrality"
        egv = nx.eigenvector_centrality_numpy(graph)
        centralities["eigenvector"] = egv
        elem_types["eigenvector"] = "float"
        cent_keys.append("eigenvector")
    if opts.flg_communicability:
        print "communicability centrality"
        cmn = nx.communicability_centrality(graph)
        centralities["communicability"] = cmn
        cent_keys.append("communicability")
    if opts.flg_load:
        print "load centrality"
        if opts.key_weight: print "weight: " + opts.key_weight
        lod = nx.load_centrality(graph, weight=opts.key_weight)
        centralities["load"] = lod
        cent_keys.append("load")

    for ck in cent_keys:
        elem_types[ck] = "float"

    ## set types for the case of element list is not given
    if len(elems.keys()) == 0:
        keys = []
        for name, d in centralities.items(): keys.extend(d.keys())
        keys = set(keys)
        for key in keys:
            elems[key] = {}
        elem_types["node_id"] = headers[0][1]

    for centrality_type, d in centralities.items():
        for elem_id, value in d.items():
            elems[elem_id][centrality_type] = value

    ## remove elements with N/A value
    new_elems = {}
    for elem_id, e in elems.items():
        flg_na = False
        for key, val_type in elem_types.items():
            if not key in e:
                flg_na = True
        if not flg_na:
            new_elems[elem_id] = e
    
    binrel = kkbinrel.BinaryRelation(symmetry=True, elem_key=opts.key_element)
    binrel.elements = new_elems
    binrel.push_type_defs(elem_types)
    print "Elems: " + str(len(binrel.elements.keys()))


    ##output_nodes(opts.fn_out, node_cent, keys, opts.flg_out_header)
    kktable.TableWriter(opts.fn_out).write_binrel_elements(binrel)

    ## grouping 
    if opts.key_group:
        key_group = opts.key_group
        g_elems, g_elem_types = kktable.TableReader(opts.fn_groups).read_binrel_elements(key_group)
    
        for ck in cent_keys:
            g_elem_types[ck] = "float"

        for cent_key in cent_keys:
            for elem_id, e in binrel.elements.items():
                if not cent_key in g_elems[e[key_group]] or\
                        g_elems[e[key_group]][cent_key] < e[cent_key]:
                    g_elems[e[key_group]][cent_key] = e[cent_key]
        g_binrel = kkbinrel.BinaryRelation(symmetry=True, elem_key="node_id")
        g_binrel.elements = g_elems
        g_binrel.push_type_defs(g_elem_types)
        kktable.TableWriter(opts.fn_out_group).write_binrel_elements(g_binrel)    
                      
    return 

def output_nodes(fn_out, nodes, keys, na_string="na", flg_out_header=False):
    f = open(fn_out,"w")
    if flg_out_header: f.write("\t".join(keys)+"\n")
    for node_id, d in nodes.items():
        line = str(node_id)
        for key in keys:
            if key in d:
                line += "\t" + str(d[key])
            else:
                line += "\t" + na_string
        line += "\n"
        f.write(line)
    f.close()

def reshape_dict(dict_of_dict):
    keys = []
    for name, d in dict_of_dict.items(): keys.extend(d.keys())
    keys = set(keys)
    
    nodes = {}
    for key in keys:
        nodes[key] = {}
        for name, d in dict_of_dict.items():
            if key in d:
                nodes[key][name] = d[key]
    return nodes

def rcp(val):
    return 1.0 / val

def add_field(graph, key_abs, func, pref):
    for node1, edges1 in graph.edge.items():
        for node2, edge in edges1.items():
            for key in key_abs:
                edge[pref+key] = func(edge[key])
    return graph

def read_headers(fn_table):
    f = open(fn_table)
    header_line = re.compile("\s+").split(f.readline().strip())
    print header_line
    f.close()
    headers = []
    for head in header_line:
        terms = re.compile("\.").split(head)
        if len(terms) > 1:
            field_name = ".".join(terms[0:-1])
            field_type = str
            if terms[-1] == "int":
                field_type = int
            elif terms[-1] == "float":
                field_type = float
            headers.append((field_name, field_type))
    return tuple(headers), header_line[0][0]


if __name__ == "__main__":
    _main()
    print "Completed."









