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
import kktrajtrans as kktt

def define_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_traj',
                 help="file name for the trajectory")
    p.add_option('--init', dest='fn_init',
                 help="file name for the initial structure")
    #p.add_option('--dt-traj', dest='dt_traj',
    #             type="float",
    #             help="timestep of frames of trajectory file [ps]")
    #p.add_option('--start-time', dest='start_time',
    #             default=0.0,
    #             type="float",
    #             help="time of the first frame")

    p.add_option('--select', dest='select',
                 #default = "name CA",
                 help="Selection of atoms to be aligned")
    p.add_option('--frame-begin', dest='frame_begin',
                 type="int", default=0,
                 help="Frame number that specifies the first frame to be considered")
    p.add_option('--frame-end', dest='frame_end',
                 type="int", default=-1,
                 help="Frame number that specifies the last frame to be considered")

    #### JOB CONTROL
    p.add_option('--n-div', dest='n_div',
                 default = 1,
                 type="int",
                 help="Number of division of atoms. The actual number of tasks will be (n-div * n-div - n-div) / 2")
    p.add_option('--task-id', dest='task_id',
                 default = 1,
                 type="int",
                 help="Task-id of this task")
    p.add_option('-o', dest='fn_output',
                 help="output file name")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def select_atoms_bynum(univ,num):
    sel_str = "bynum "+" or bynum ".join(num)
    atoms = univ.selectAtoms(sel_str)
    return atoms

def cal_distances_pairs(kk_traj_trans, fn_out,
                        atoms, pair_begin, pair_end,
                        trajreader, frame_begin, frame_end):
    n_atoms = len(atoms)
    f_out = open(fn_out,'w')
    headers = ["atom_id1.int", "atom_id2.int", "ave.float", "sd.float"]
    #print "\t".join(headers)+"\n"
    cur_pair = -1
    for atom1 in atoms:
        atom_id1 = atom1.number
        crd1 = trajreader.read_atom_crd(atom_id1)
        for atom2 in atoms:
            atom_id2 = atom2.number
            if atom_id1 >= atom_id2: continue
            cur_pair += 1
            if pair_begin > cur_pair: continue
            if pair_end <= cur_pair: break

            crd2 = trajreader.read_atom_crd(atom_id2)
            crd_d = crd1-crd2
            ##crd_d[0]*crd_d[0] + crd_d[1]*crd_d[1] + crd_d[2]*crd_d[2]
            dist = np.sqrt(np.sum(crd_d * crd_d, axis=0))
            if frame_begin < 0: frame_begin = 0
            if frame_end < 0:
                dist = dist[frame_begin:]
            else:
                dist = dist[frame_begin:frame_end]

            ave = np.mean(dist) * 10.0
            std = np.std(dist)  * 10.0
            min_d = np.min(dist) * 10.0
            max_d = np.max(dist) * 10.0

            line = str(atom_id1) + "\t" + str(atom_id2)
            line += "\t" + str(ave)
            line += "\t" + str(std)
            line += "\t" + str(min_d)
            line += "\t" + str(max_d)
            line += "\n"
            f_out.write(line)
    f_out.close()

def _main():
    opts, args = define_options() 
    trajreader = kktt.TrajTReader(opts.fn_traj)
    trajreader.open()
    n_atoms, n_frames = trajreader.read_header()

    fn_ref = opts.fn_init
    univ = Universe(fn_ref)
    sel_atoms = univ.selectAtoms(opts.atom_select)

    print("Task:")
    n_sel_atoms = len(sel_atoms)
    n_pairs = ((n_sel_atoms * n_sel_atoms) - n_sel_atoms)/2
    task_id = opts.task_id - 1

    n_pairs_per_task = (n_pairs + ( opts.n_div - 1)) / opts.n_div
    print(str(task_id+1) + " / " + str(opts.n_div))


    #print "<<<ATOM>>>"
    #print "num.int atom_id.int atom_name.str res_id.int res_num.int res_name.str"
    #for i, atom in enumerate(sel_atoms):
    #    line = str(i) + " " + str(atom.number) +" " + str(atom.name) + " "
    #    line += str(atom.resid) + " " + str(atom.resnum) + " " + str(atom.resname)
    #    print line
    print ""

    pair_begin = n_pairs_per_task * task_id
    pair_end = pair_begin + n_pairs_per_task

    print "cal_distances"
    cal_distances_pairs(opts.fn_traj, opts.fn_output, 
    sel_atoms, pair_begin, pair_end,
                        trajreader, opts.frame_begin, opts.frame_end)
    trajreader.close()
    print "Terminated"
    
if __name__ == "__main__":
    _main()

class Test(object):
    ## test code
    ##  four atoms
    ##  1 ns 
    def __init__():
        self.fn_traj = "test_traj.dat"
        self.fn_pass = "test_pass.txt"
        self.fn_state = "test_state.txt"
        self.traj = None
        return 
    def test_gen_traj(self):
        return
    def test_gen_pass(self):
        return
    def test_gen_state(self):
        return
    
