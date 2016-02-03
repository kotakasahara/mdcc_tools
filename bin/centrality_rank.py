#!/usr/bin/python2.7

from optparse import OptionParser
import kkkit
import kktext
import kkdefine
import kktable
import kkbinrel
import re
import sys

def define_options():
    p = OptionParser()
    
    p.add_option('--i-node', dest='fn_i_nodes',
                 help="file name for nodes")
    p.add_option('--i-edge', dest='fn_i_edges',
                 help="file name for edges")
    p.add_option('--o-node', dest='fn_o_nodes',
                 help="file name for nodes")
    p.add_option('--o-edge', dest='fn_o_edges',
                 help="file name for edges")

    p.add_option('--key-elem', dest='key_element',
                 default="node_id",
                 help="key of elements")
    p.add_option('--key-src', dest='key_src',
                 default="node_id",
                 help="key of edge source column")
    p.add_option('--key-dst', dest='key_dst',
                 default="node_id",
                 help="key of edge destination column")

    p.add_option('--desc', dest='flg_desc',
                 action="store_true",
                 help="Descending order")
    p.add_option('--col-value', dest='col_value',
                 help="Column ID for ranking measure")
    p.add_option('--col-category', dest='col_category',
                 help="Column ID for ranking category")
    p.add_option('--category', dest='category',
                 action="append",
                 help="Category; values of col_category concatenating with the  symbol '|' ")

    p.add_option('--rank-col', dest='rank_column',
                 default="rank",
                 help="Name of the newly added column for ranking information")


    p.add_option('--header', dest='flg_header',
                 action="store_true",
                 help="Header")
        
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def rank_category(binrel, cate,
                  col_category, col_value,
                  flg_desc):
    values = {}
    #print binrel.elements.keys()[1:30]
    
    for elem_key, elem in binrel.elements.items():
        if elem[col_category] in cate:
            values[elem_key] = elem[col_value]
            #print elem

    values_st = sorted(values.items(), key=lambda x:x[1], reverse=flg_desc)

    rank = {}
    for i, kv in enumerate(values_st):
        rank[kv[0]] = i

    return rank

def _main():
    opts, args = define_options() 

    trm = kktable.TableReadManager(opts.fn_i_nodes, 
                                   opts.fn_i_edges,
                                   opts.key_element,
                                   opts.key_src, opts.key_dst)
    tbl = trm.read_tables(False)

    rank_column = opts.rank_column
    rank_column_1 = opts.rank_column+"1"
    rank_column_2 = opts.rank_column+"2"

    tbl.push_type_def(rank_column, "int")
    tbl.push_type_def(rank_column_1, "int")
    tbl.push_type_def(rank_column_2, "int")

    ## ranking each node
    for cate in opts.category:
        set_cate = set(cate.split(","))
        print set_cate
        rank = rank_category(tbl, set_cate,
                             opts.col_category, opts.col_value,
                             opts.flg_desc)
        print len(rank.keys())
        print rank.values()

        for elem_key, rank in rank.items():
            tbl.elements[elem_key][rank_column] = rank
            
    ## adding ranking information to edges
    for edge_keys, edge_info in tbl.relations.items():

        if not rank_column in tbl.elements[edge_keys[0]]:
            print "ERROR"
            print edge_keys
            print tbl.elements[edge_keys[0]]
            #sys.exit(1)
            edge_info[rank_column_1] = -1
        else:
            edge_info[rank_column_1] = tbl.elements[edge_keys[0]][rank_column]

        if not rank_column in tbl.elements[edge_keys[1]]:
            print "ERROR"
            print edge_keys
            print tbl.elements[edge_keys[1]]
            #sys.exit(1)
            edge_info[rank_column_2] = -1
        else:
            edge_info[rank_column_2] = tbl.elements[edge_keys[1]][rank_column]
    
    kktable.TableWriter(opts.fn_o_nodes).write_binrel_elements(tbl)

    kktable.TableWriter(opts.fn_o_edges).write_binrel_relations(tbl)

if __name__ == "__main__":
    _main()
