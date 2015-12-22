#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import numpy as np
import scipy.stats
from optparse import OptionParser
from MDAnalysis import Universe
from MDAnalysis.analysis import distances
import re
import MDAnalysis
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False
import sys
import os

def define_options():
    p = OptionParser()
    p.add_option('-i', '--i-pdb', dest='fn_pdb',
                 help="file name for structure file in .pdb format")
    p.add_option('-s', '--select', dest='select',
                 default="*",
                 help="selection string in MDAnalysis syntax")
    p.add_option('-o', '--o-dist', dest='fn_dist',
                 help="file name for output distance map")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = define_options() 
    univ = Universe(opts.fn_pdb)
    sel_atoms = univ.selectAtoms(opts.select)
    print "n_atoms: " + str(len(sel_atoms))
    f_out = open(opts.fn_dist, "w")
    for a1 in sel_atoms:
        for a2 in sel_atoms:
            dist = np.linalg.norm(a1.pos - a2.pos)
            line = str(a1.number) + "\t" + str(a2.number) + "\t"
            line += str(dist) + "\n"
            f_out.write(line)
    f_out.close()
            
if __name__ == "__main__":
    _main()
