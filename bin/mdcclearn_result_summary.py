#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd
from optparse import OptionParser
import re
import os

def define_options():
    p = OptionParser()
    
    p.add_option('--dir-mdcclearn', dest='dir_mdcclearn',
                 help="directory containing mdcclearn outputs.")
    p.add_option('--pref-mdcclearn', dest='pref_mdcclearn',
                 help="Prefix for mdcclearn results outputs.")
    
    p.add_option('-o', dest='fn_out',
                 help="output file name")
    p.add_option('--min-pi', dest='min_pi',
                 type="float", default=0.01,
                 help="minimum threshold of pi")
    p.add_option('--dim', dest='dimension',
                 type="int", default=3,
                 help="Dimension")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def get_header_string(dim):
    header = ["gc_id.int", "element_id.int", "pi.float"]
    for i in range(dim):
        header.append("mu"+str(i+1)+".float")
    for i in range(dim):
        for j in range(dim):
            header.append("sigma"+str(i+1)+str(j+1)+".float")
    return "\t".join(header)

def _main():
    opts, args = define_options() 
    f_out = open(opts.fn_out,"w")
    f_out.write(get_header_string(opts.dimension) + "\n")
    filelist = {}
    for fn1 in os.listdir(opts.dir_mdcclearn):
        m = re.compile(opts.pref_mdcclearn + "(\d+)").match(fn1)
        if m:
            type_id = m.group(1)
            filelist[int(type_id)] = os.path.join(opts.dir_mdcclearn, fn1)

    gc_id = 0
    for i, type_id in enumerate(sorted(filelist.keys())):
        #print str(type_id) + " : " + str(i) + " / " + str(len(filelist.keys()))
        f = open(filelist[type_id])
        for line in f:
            terms = re.compile("\s+").split(line.strip())
            if float(terms[1]) < opts.min_pi: continue
            new_line = [str(gc_id), str(type_id)]
            new_line.extend(terms[1:])
            f_out.write("\t".join(new_line)+"\n")
            gc_id += 1
        f.close()
    f_out.close()
    print "finished."
    return

if __name__ == "__main__":
    _main()
