#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import numpy as np
from optparse import OptionParser
import MDAnalysis as mda
import re
import struct as st
import sys


def define_options():
    p = OptionParser()
    p.add_option('--i-pdb', dest='fn_pdb',
                 help="The structure file in .pdb format")
    p.add_option('-i', dest='fn_trr',
                 action="append",
                 help="Trajectory file / .trr")
    p.add_option('--i-list', dest='fn_list',
                 help="file list")
    p.add_option('-o', dest='fn_out',
                 help="Output trajectory file")
    p.add_option('-f', dest='out_format',
                 type="choice",
                 choices=["ncdf", "trr", "dcd"],
                 help="Output trajectory format")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args


def check_options(opts):
    if opts.fn_list and opts.fn_trr:
        sys.stderr.write("Do not specify both -i and --i-list\n")
        sys.exit(1)
    if not opts.fn_list and not opts.fn_trr:
        sys.stderr.write("One of -i or --i-list are required.\n")
        sys.exit(1)
    return 

def read_file_list(fn_list):
    files = []
    f = open(fn_list)
    for line in f:
        files.append(line.strip().split()[0])
    return files

def check_file_format(files):
    ext = set()
    for fn in files:
        ext.add(fn[-3:])
    if len(ext) > 1:
        sys.stderr.write("The file format of the trajectory files should be consistent.\n")
        sys.stderr.write(", ".join(ext)+"\n")
        sys.exit(1)
    return list(ext)[0]

def convert(fn_out, fn_pdb, files, trj_format, out_format):
    u = mda.Universe(fn_pdb, files)
    writer = None
    if out_format == "ncdf":
        writer = mda.coordinates.TRJ.NCDFWriter(fn_out, len(u.atoms))
    elif out_format == "trr":
        writer = mda.coordinates.TRR.TRRWriter(fn_out, len(u.atoms))
    elif out_format == "dcd":
        writer = mda.coordinates.DCD.DCDWriter(fn_out, len(u.atoms))
    for ts in u.trajectory:
        writer.write_next_timestep(ts)
    return 

def _main():
    opts, args = define_options() 
    check_options(opts)
    files = []
    if opts.fn_list:
        files = read_file_list(opts.fn_list)
    else:
        files = opts.fn_trr
        
    trj_format = check_file_format(files)
    convert(opts.fn_out, opts.fn_pdb, files, trj_format, opts.out_format)

    return 


if __name__ == "__main__":
    _main()

