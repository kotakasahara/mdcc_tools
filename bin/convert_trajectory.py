#!/usr/bin/python2.7
#$ -S /usr/bin/python2.7
#$ -cwd

import numpy as np
from optparse import OptionParser
from MDAnalysis import Universe
import re
import struct as st
import sys

def define_options():
    p = OptionParser()
    p.add_option('-i', dest='fn_trr',
                 action="append",
                 help="Trajectory file / .trr")
    p.add_option('--i-list', dest='fn_list',
                 help="file list")
    p.add_option('-o', dest='fn_out',
                 help="Output trajectory file")
    p.add_option('--i-format', dest='in_format',
                 type="choice",
                 choices=["trr", "tsv"],
                 default="trr",
                 help="Input file format, tsv or trr")
    p.add_option('--ignore-col', dest='ignore_col',
                 type="int", default=0,
                 help="The number of ignored columns in .tsv files.")
    p.add_option('--ignore-row', dest='ignore_row',
                 type="int", default=0,
                 help="The number of ignored rows in .tsv files.")
    p.add_option('--double', dest='flg_double',
                 action="store_true",
                 help="Flag for double precision")
    p.add_option('--dim', dest='dim',
                 type="int", default=3,
                 help="Dimension")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def check_options(opts):
    if opts.in_format == ".trr":
        if not opts.fn_trr or len(opts.fn_trr) > 1 or fn_list:
            print "Only one .trr file must be specified by -i option"
            sys.exit(1)
    return

class Traj(object):
    def __init__(self, fn_list, fn_out, i_fmt, dim, ignore_row, ignore_col):
        self.fn_list = fn_list
        #self.f = open(fn,"rb")
        self.fo = open(fn_out,"wb")
        self.real_code = "f"
        self.size_real = 4
        self.endian = ">" #big endian
        self.dim = dim
        self.in_format = i_fmt
        self.ignore_row = ignore_row
        self.ignore_col = ignore_col
    def open_trr(self):
        self.f = open(self.fn_list[0],"rb")
    def close(self):
        self.f.close()
        self.fo.close()
    def set_endian(self):
        #Checking endian
        self.f.seek(0)
        buf = self.f.read(4)
        ##print buf
        ##print len(buf)
        ##print str(st.unpack("<i",buf)) + ":"+str(st.unpack(">i",buf))
        if st.unpack("<i",buf)[0] == 1993:
            self.endian = "<"
            print("Little endian")
        elif st.unpack(">i",buf)[0] != 1993:
            sys.stderr.write("Error: The first 4 byte is not the integer 1993")
        else:
            self.endian = ">"            
            print("Big endian")
        return  self.endian
    def set_real(self, flg_double):
        if flg_double:
            self.code_real = "d"
            self.size_real = 8
        else:
            self.code_real = "f"
            self.size_real = 4
    def read_int(self):
        """
        unpack a value from binary string
        """
        read_string = self.endian + "i"
        i = st.unpack(read_string, self.f.read(4))[0]
        return i
    def read_real(self):
        """
        unpack a value from binary string
        """
        read_string = self.endian + self.code_real
        i = st.unpack(read_string, self.f.read(self.size_real))[0]
        return i
    def write_int(self, val):
        write_string = self.endian + "i"
        buf = st.pack(write_string, val)
        self.fo.write(buf)
        return buf
    def get_bin_real(self, val):
        write_string = self.endian + self.code_real
        bin = st.pack(write_string, val)
        return bin
    def write_real(self, val):
        buf = get_bin_real(val)
        self.fo.write(buf)
        return buf
    def read_information_tsv(self):
        """
        read basic information from the tsv files
        self.n_atoms
        self.n_frames
        self.atom_dict
          atom_dict[atom_id] = [file_id, elem_id_in_each_file]
        """

        ## count self.n_frames
        self.n_frames = 0
        f = open(self.fn_list[0])
        print self.fn_list[0]
        for i in range(self.ignore_row):
            f.readline()
        for line in f:
            self.n_frames += 1
        f.close()            
        print "n_frames : " + str(self.n_frames)

        # set self.n_atoms, atom_dic
        self.n_atoms = 0
        self.atom_dict = {}
        ## atom_dic
        for i_f, fn in enumerate(self.fn_list):
            f = open(fn)
            for i in range(self.ignore_row):
                f.readline()

                #for i_f, line in enumerate(f):
            line = f.readline()
            terms = line.strip().split()[self.ignore_col:]
            n_atom_f = len(terms)/self.dim
            for i_atom in range(n_atom_f):
                self.atom_dict[self.n_atoms] = (i_f, i_atom)
                self.n_atoms += 1
            f.close()
        print "n_atoms : " + str(self.n_atoms)
        print self.atom_dict
        return

    def read_information(self):
        """
        read basic information of this trajectory file
        self.box_size
        self.x_size
        self.v_size
        self.f_size
        self.n_atoms
        self.n_frames
        """
        self.f.seek(0)
        self.f.read(32)
        self.size_box = self.read_int() #32-36
        self.f.read(16)                 #36-52
        self.size_x = self.read_int()   #52-56
        self.size_v = self.read_int()   #56-60
        self.size_f = self.read_int()   #60-64
        self.n_atoms = self.read_int()  #64-68
        self.f.read(8) #68-76
        self.read_real() #time
        self.read_real() #lambda
        
        self.size_header = 84
        if self.code_real=="d": 
            self.size_header = 92
        
        self.size_frame = self.size_header + self.size_box + self.size_x + self.size_v + self.size_f
        
        ### ?????????????  2013.11.22
        ####if self.code_real=="d": self.size_frame -= 8
        ### ?????????????  2013.11.22
        
        self.f.seek(0)
        self.n_frames = -1
        buf = "dummy"
        while buf != "":
            self.n_frames += 1
            self.f.seek(self.n_frames * self.size_frame)
            buf = self.f.read(4)
        self.f.seek(0)
        print "n_frames: " + str(self.n_frames)
        print "n_atoms: " + str(self.n_atoms)
        print "size_x: " + str(self.size_x)
        print "size_v: " + str(self.size_v)
        print "size_f: " + str(self.size_f)
        print "size_frame: " + str(self.size_frame)
        return

    def read_atom_pos(self,i_atom):
        if self.in_format=="trr":
            return self.read_atom_pos_trr(i_atom)
        elif self.in_format=="tsv":
            return self.read_atom_pos_tsv(i_atom)            
    def read_atom_pos_tsv(self,i_atom):
        f = open(self.fn_list[self.atom_dict[i_atom][0]])
        dat = []
        for d in range(self.dim):
            dat.append("")
        for i in range(self.ignore_row):
            f.readline()
        for line in f:
            terms = line.strip().split()[self.ignore_col:]
            for d in range(self.dim):
                val = float(terms[self.atom_dict[i_atom][1]*self.dim+d])
                dat[d] += st.pack(self.endian+"f", val)
        f.close()
        return dat
    def read_atom_pos_trr(self,i_atom):
        crd_x = ""
        crd_y = ""
        crd_z = ""
        if self.size_x > 0:
            pos_init = self.size_header + self.size_box + i_atom * self.size_real * 3
            self.f.seek(pos_init)
            for i_frame in xrange(self.n_frames):
                #pos = pos_init + i_frame * self.size_frame
                #print str(i_atom) + " " + str(pos) + "/" + str(self.size_frame * self.n_frames)
                #self.f.seek(pos)
                tmp_x = self.f.read(self.size_real)
                tmpf_x = st.unpack(self.endian+self.real_code, tmp_x)[0]*10.0
                crd_x += st.pack(self.endian+"f",tmpf_x)

                tmp_y = self.f.read(self.size_real)
                tmpf_y = st.unpack(self.endian+self.real_code, tmp_y)[0]*10.0
                crd_y += st.pack(self.endian+"f",tmpf_y)

                tmp_z = self.f.read(self.size_real)
                tmpf_z = st.unpack(self.endian+self.real_code, tmp_z)[0]*10.0
                crd_z += st.pack(self.endian+"f",tmpf_z)

                self.f.seek(self.size_frame-self.size_real*3, 1)
        self.f.seek(0)
        return crd_x, crd_y, crd_z
    def write_header(self):
        #0-3: magic_number 1993
        #4-7: 8:double or 4:float
        #8-11: n_atoms
        #12-15: n_frames
        #16-19: dimension
        self.write_int(1993)
        self.write_int(self.size_real)
        self.write_int(self.n_atoms)
        self.write_int(self.n_frames)
        self.write_int(self.dim)
        return
    def write_atom(self,crd):
        assert(len(crd) == self.dim)
        for i_crd in crd:
            self.fo.write(i_crd)
        return 

def read_fn_list(fn_fn_list):
    fn = []
    f = open(fn_fn_list)
    for line in f:
        terms = line.strip().split()
        fn.append(terms[0])
    return fn

def _main():
    """
    Generating mDCC trajectory file from the Gromacs .trr file or tsv
    Output a new binary file as below:
    28 bytes for File header
      0-3    : the magic number "1993"
      4-7    : size_real (4 or 8)
      8-11   : number of elements
      12-15  : number of frames
      16-19  : dimension
    Following parts records the coordinates
      xxxx...yyyy...zzzz...xxxx...yyyy...
    """
    
    opts, args = define_options() 
    check_options(opts)
    read_pref = ">" # big endian

    print("read trajectory")
    if opts.fn_trr:
        fn_list = opts.fn_trr
    elif opts.fn_list:
        fn_list = read_fn_list(opts.fn_list)

    traj = Traj(fn_list, opts.fn_out, opts.in_format,
                opts.dim,
                opts.ignore_row, opts.ignore_col)

    if opts.in_format == "trr":
        traj.open_trr()
        traj.set_endian()
        traj.set_real(opts.flg_double)
        traj.read_information()
    elif opts.in_format == "tsv":
        traj.read_information_tsv()
        
    traj.write_header()
    
    for i_atom in range(traj.n_atoms):
        #print "Atom : " + str(i_atom)
        xyz = traj.read_atom_pos(i_atom)
        traj.write_atom(xyz)

    if opts.in_format == "trr":
        traj.close()
    return 

if __name__ == "__main__":
    _main()
