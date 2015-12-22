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
                 help="Trajectory file / .trr")
    p.add_option('-o', dest='fn_out',
                 help="Output trajectory file")
    p.add_option('--double', dest='flg_double',
                 action="store_true",
                 help="Flag for double precision")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args


class Traj(object):
    def __init__(self, fn, fn_out):
        self.f = open(fn,"rb")
        self.fo = open(fn_out,"wb")
        self.real_code = "f"
        self.size_real = 4
        self.endian = ">" #big endian
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
                crd_x += self.f.read(self.size_real)
                crd_y += self.f.read(self.size_real)
                crd_z += self.f.read(self.size_real)
                self.f.seek(self.size_frame-self.size_real*3, 1)
        self.f.seek(0)
        return crd_x, crd_y, crd_z
    def write_header(self):
        #0-4: magic_number 1993
        #4-8: 8:double or 4:float
        #8-12: n_atoms
        #12-16: n_frames
        self.write_int(1993)
        self.write_int(self.size_real)
        self.write_int(self.n_atoms)
        self.write_int(self.n_frames)
        return
    def write_atom(self,crd_x,crd_y,crd_z):
        self.fo.write(crd_x)
        self.fo.write(crd_y)
        self.fo.write(crd_z)
        return 

def _main():
    """
    Transpositioning Gromacs .trr file.
    Output a new binary file as below:
    24   : File header
      0-4    : the magic number "1993"
      4-8    : size_box
      8-12   : number of frames
      12-16  : size_x (number of atoms)
      16-20  : size_v (number of atoms)
      20-24  : size_f (number of atoms)
    : Frame
    0-4: step
    4-8: time
    size_box: 
    number of frames: x coordinates of the first atom
    number of frames: y coordinates of the first atom
    number of frames: z coordinates of the first atom
    number of frames: x coordinates of the second atom
    ...
    number of frames: z coordinates of the last atom 
    number of frames: x velocities of the first atom   
    ...
    """
    
    opts, args = define_options() 
    read_pref = ">" # big endian


    print("read trajectory")
    traj = Traj(opts.fn_trr, opts.fn_out)
    traj.set_endian()
    traj.set_real(opts.flg_double)
    traj.read_information()

    traj.write_header()
    
    for i_atom in range(traj.n_atoms):
        #print "Atom : " + str(i_atom)
        x,y,z = traj.read_atom_pos(i_atom)
        traj.write_atom(x,y,z)
        
    traj.close()
    return 
if __name__ == "__main__":
    _main()
