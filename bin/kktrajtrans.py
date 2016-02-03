#!/usr/bin/python2.7

import struct as st
import kkkit
import numpy as np

class TrajTReader(kkkit.FileBI):
    def __init__(self, fn):
        super(TrajTReader, self).__init__(fn)
        self.size_real = 8
        self.fc_read_real = None
        self.n_atoms = None
        self.n_frames = None
        self.dim = None
        self.size_header = 20
    def open(self):
        return super(TrajTReader, self).open()
    def read_header(self):
        self.f.seek(0)
        buf = self.f.read(4)
        if st.unpack("<i",buf)[0] == 1993:
            self.endian = "<" #little endian
        elif st.unpack(">i",buf)[0] == 1993:
            self.endian = ">" #big endian
        else:
            sys.stderr.write("Error: The first 4 byte is not the integer 1993")
            self.endian = ""            
        
        ## size of real values [4 or 8 ]
        self.size_real = self.read_int()[0]
        if self.size_real == 4:
            self.fc_read_real = self.read_float
        elif self.size_real == 8:
            self.fc_read_real = self.read_double
        ## number of atoms
        self.n_atoms = self.read_int()[0]
        ## number of frames
        self.n_frames = self.read_int()[0]
        ## dimension
        self.dim = self.read_int()[0]
        return self.n_atoms, self.n_frames
    def read_atom_crd(self, atom_id, read_from=0, read_to=-1 ):
        ## begin_from ... number of frame specifying the point of begining to read
        if not 1 <= read_to <= self.n_frames:
            read_to = self.n_frames
        n_read_frames = read_to - read_from 

        crd = np.zeros((self.dim, self.n_frames), dtype="float32")
        
        for d in range(self.dim):
            self.f.seek(self.size_header + (atom_id * self.dim + d) * self.n_frames * self.size_real + read_from * self.size_real)
            crd[d] = np.array(self.fc_read_real(n_read_frames))

        return crd
    
        
