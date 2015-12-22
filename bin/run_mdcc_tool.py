#!/usr/bin/python2.7

from optparse import OptionParser
import re
from MDAnalysis import Universe

def define_options():
    p = OptionParser()
    p.add_option('--pdb', dest='fn_pdb',
                 help="Structure file")
    p.add_option('--select', dest='selection',
                 help="Atom selection string for MDAnalysis")
    p.add_option('--n-div', dest='n_div',
                 type="int",
                 help="Total number of jobs")
    p.add_option('--task-id', dest='task_id',
                 type="int",
                 help="Task ID of this task (begin with 1)")
    p.add_option('--mdcc-conf', dest='fn_mdcc_conf',
                 help="Template file for the configuration")
    p.add_option('-o', dest='pref_out',
                 help="output file prefix")
    p.add_option('--fn-sh', dest='fn_sh',
                 default="run_mdcc.bash",
                 help="intermediate sh filename")
    p.add_option('--mdcc-bin', dest='mdcc_bin',
                 help="Path to the binary file of a mdcc program (mdcc_learn or mdcc_assign).")
    p.add_option('--mode', dest="mode",
                 type="choice",
                 choices=["mdcclearn","mdccassign"],
                 default="mdcclearn",
                 help="option keyword for data table file name")
    p.add_option('--table-bin' , dest="fn_table",
                 help="binary tablefile")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def read_configs(fn_conf):
    f = open(fn_conf)
    conf = ""
    for line in f:
        conf += line.strip() + " "
    f.close()
    return conf

def _main():
    """
    
    """
    opts, args = define_options() 
    mdcclearn_conf = read_configs(opts.fn_mdcc_conf)
    f_bash = open(opts.fn_sh,'w')

    f_bash.write("#!/bin/bash\n")
    f_bash.write("#$ -S /bin/bash\n")
    f_bash.write("#$ -cwd\n")
    #f_bash.write("source "+SHELLRC + "\n")
    f_bash.write("\n")

    univ = Universe(opts.fn_pdb)
    atoms = univ.selectAtoms(opts.selection)
    n_atoms_task = (len(atoms) + (opts.n_div - 1)) / opts.n_div
    task_id = opts.task_id - 1
    atom_begin = n_atoms_task * task_id
    atom_end = min(len(atoms), atom_begin + n_atoms_task)

    print "The number of selected atoms : %d"%len(atoms)
    print "The number of atoms processed in this job : %d (%d - %d)"%(n_atoms_task,atom_begin, atom_end)
    
    for i in range(atom_begin, atom_end):
        atom_id = atoms[i].number
        
        conf = mdcclearn_conf.replace("#{COLUMN}",str(atom_id))

        cmd = opts.mdcc_bin + " " + conf 

        f_bash.write(cmd + '\n')
    f_bash.close()

def dummy_task(f_bash, mode, mdcc_bin, conf):
    fn_dummy_table = "dummy_table.dat"
    fn_dummy_out = "dummy_out.txt"
    f_dum_table = open(fn_dummy_table,"w")
    for i in range(0,10):
        f_dum_table.write("\t".join([str(x) for x in range(0,100)]))
    f_dum_table.write("\n")
    f_dum_table.close()
    cmd = ""
    if mode=="mdccassign":
        table_kw = "-fn-interactions"
        conf += " " + table_kw + " " + fn_dummy_table
        conf += " -fn-result " + fn_dummy_out
        conf += " -gmm-type " + "0"
        cmd = mdcc_bin + " " + conf 
    elif mode=="mdcclearn":
        table_kw = "-fn-data-table"
        conf += " " + table_kw + " " + fn_dummy_table
        cmd = mdcc_bin + " " + conf + " > " + fn_dummy_out
    f_bash.write("\n" + cmd + '\n')    
    return

if __name__ == "__main__":
    _main()
