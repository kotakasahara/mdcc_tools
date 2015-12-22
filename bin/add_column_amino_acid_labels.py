#!/usr/bin/python2.7

from optparse import OptionParser
import kkkit
import kktext
import kkdefine
import re

def define_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_inp',
                 help="file name for input table")
    p.add_option('-o', dest='fn_out',
                 help="file name for output table")
    p.add_option('--col-res', dest='col_res',
                 type="int",
                 help="Column ID for residues name")
    p.add_option('--header', dest='flg_header',
                 action="store_true",
                 help="Header")
        
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = define_options() 
    if not opts.fn_inp:
        return
    f = open(opts.fn_inp)
    fo = open(opts.fn_out,"w")
    if opts.flg_header:
        line = f.readline()
        terms = re.compile("\s+").split(line.strip())
        terms.append("res_char.str")
        fo.write("\t".join(terms)+"\n")
    for line in f:
        terms = re.compile("\s+").split(line.strip())        
        ch = "X"
        if terms[opts.col_res][0:3] in kkdefine.KKDEF.AA_3_1:
            ch = kkdefine.KKDEF.AA_3_1[terms[opts.col_res][0:3]]
        elif terms[opts.col_res][0:3] in kkdefine.KKDEF.NA_3_1:
            ch = kkdefine.KKDEF.NA_3_1[terms[opts.col_res][0:3]]
        terms.append(ch)
        fo.write("\t".join(terms)+"\n")            
    f.close()
    fo.close()
    
if __name__ == "__main__":
    _main()
