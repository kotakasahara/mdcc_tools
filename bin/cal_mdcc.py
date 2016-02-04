#$ -S /usr/bin/python2.7
#$ -cwd

from optparse import OptionParser
from MDAnalysis import Universe
import sys
import numpy as np
import os
import re
import struct
import kktrajtrans
import local_block_bootstrap as lbb

DEBUG = True

EPS6 = 1e-8

def define_options():
    """
    Three files are needed.
    [--dir-res]_[res_num]/[--pref-assign]_[atom_id].txt
    [--dir-res]_[res_num]/[--pref-crd]_[atom_id].txt
    """

    p = OptionParser()

    p.add_option('--fn-ref', dest='fn_ref',
                 help="file name for ref file")
    p.add_option('--skip', dest='skip',
                 default=1,
                 type="int",
                 help="number of steps to be skipped")
    #p.add_option('--center-mode', dest='center_mode',
    #             default="default",
    #             help="definition of the center of fluctuation")
    p.add_option('--gaussian', dest='fn_gaussian',
                 help="file name for definition of gaussian mixtures")

    p.add_option('--assign-binary', dest='flg_assign_bin',
                 action="store_true",
                 help="The flag indicating that the assignment files were written in binary format")
    p.add_option('--pref-assign', dest='pref_assign',
                 help="prefix of file name for trajectory of assignments to GMM.")
    p.add_option('--suff-assign', dest='suff_assign',
                 default="",
                 help="suffix of file name for trajectory of assignments to GMM.")

    p.add_option('--fn-crd-bin', dest='fn_crd_bin',
                 help="file name for trajectory of coordinates.")
    p.add_option('--pref-crd', dest='pref_crd',
                 help="prefix of file name for trajectory of coordinates.")
    p.add_option('--suff-crd', dest='suff_crd',
                 default="",
                 help="suffix of file name for trajectory of coordinates.")
    p.add_option('--crd-tsv-skip-header', dest='crd_tsv_skip_header',
                 default=0,
                 help="The number of lines to be skipped in the tsv file.")
    p.add_option('--crd-tsv-columns', dest='crd_tsv_columns',
                 action="append", type="int",
                 help="The columns-IDs defines the data in the tsv file.")

    p.add_option('--timestep', dest='timestep',
                 type="float",
                 help="timestep of each frame in picosec")
    p.add_option('--start-time', dest='start_time',
                 default=0.0,
                 type="float",
                 help="time of the first frame")
    p.add_option('-b', '--frame-begin', dest='range_time_begin',
                 type="int", default=0,
                 help="The time begining to be taken into account. The frame of just this time is included, i.e., range is defined as [-b,-e). ")
    p.add_option('-e', '--frame-end', dest='range_time_end',
                 type="int", default=-1,
                 help="The time ending to be taken into account. The frame of just this time is not included, i.e., range is defined as [-b,-e). ")

    ## considering only contacting atoms
    #p.add_option('--i-dist', dest='fn_dist_summary',
    #             help="file name for dist summary")    
    #p.add_option('--max-dist', dest='max_dist',
    #             type="float", default=5.0,
    #             help="maximum value of average interatomic distance")    

    p.add_option('--select', dest='str_select',
                 action="append",
                 help="MDAnalysis selection string for atom set 1")    
    p.add_option('--select-id', dest='str_select_id',
                 help="Atom id to be analyzed. ex) 0-5,7,10-12")

    p.add_option('--o-dcc', dest="fn_o_dcc",
                 help="")
    p.add_option('--o-mdcc', dest="fn_o_mdcc",
                 help="")
    p.add_option('--o-rmsf', dest="fn_o_rmsf",
                 help="")

    p.add_option('--min-corr', dest="min_corr",
                 type="float",default=0.0,
                 help="")
    p.add_option('--min-pi', dest="min_pi",
                 type="float",default=0.0,
                 help="")

    ## Local Block Bootstrap
    #p.add_option('--lbb-block-size', dest="lbb_block_size",
    #             type="int",
    #             help="Block size for local block bootstrap")
    #p.add_option('--lbb-b', dest="lbb_b",
    #             type="float",
    #             help="B value for local block bootstrap")
    #p.add_option('--lbb-repeat', dest="lbb_repeat",
    #             type="int",
    #             help="number of bootstrap samples for local block bootstrap")

    ## paralell
    p.add_option('--n-div', dest="n_div_job",
                 type="int",default=1,
                 help="Number of tasks")
    p.add_option('--task-id', dest="task_id",
                 type="int",default=-1,
                 help="Task id")

    ## option
    p.add_option('--coef-mode', dest="coef_mode",
                 type="choice",
                 choices=["product_weight","pkpl", "pkpk"],
                 default="pkpl",
                 help="Weighting coefficient for each Gaussian pair incorporates the product term")
    ## product_weight
    ## pkpl
    ## pkpk
    ## sum

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"

    return opts, args

def print_opt_info(opts):
    if not opts.pref_assign:
        print "The conventional Dynamic Cross Correlation"
    else:
        print "The multi-modal Dinamic Cross Correlation"

    if not opts.fn_crd_bin:
        if not opts.pref_crd or \
                not opts.crd_tsv_columns:
            print "When the binary trajectory data is not specified by --fn-crd-bin optsion,"
            print "the .tsv files must be specified by --pref-crd and --suff-crd options."
            print "The names of .tsv files must include 0-origin element (atom) ID."
            print "ex) crd_0.tsv, crd_1.tsv, ..."
            print "    --pref-crd crd_  --suff-crd .tsv"
            print ""
            print "In addition, 0-origin column IDs in the .tsv file is specified by --crd-tsv-columns option,"
            print "in order to define the column of the data to be analyzed."
            sys.exit(1)

    return 

def _main():
    opts, args = define_options() 
    print_opt_info(opts)

    fn_ref = None
    univ = None

    if opts.fn_ref:
        fn_ref = opts.fn_ref
        univ = Universe(fn_ref)

    time_range = (opts.range_time_begin, opts.range_time_end)

    atom_ids = None
    #dim = None
    flg_crd_bin = opts.fn_crd_bin != None
    
    if opts.fn_crd_bin:
        trajreader = kktrajtrans.TrajTReader(opts.fn_crd_bin)
        trajreader.open()
        trajreader.read_header()
        dim = trajreader.dim
        n_frames = trajreader.n_frames
        n_atoms = trajreader.n_atoms
        trajreader.close()
    else:
        dim = len(opts.crd_tsv_columns)
        
    if opts.str_select_id:
        atom_ids = select_atom_ids(opts.str_select_id)
    elif univ:
        if len(opts.str_select)>0:
            print opts.str_select
            atom_ids = get_atom_ids_from_select(univ, opts.str_select)
        else:
            atom_ids = get_atom_ids_from_select(univ, ["name *"])
    else:
        atom_ids = range(n_atoms)

    
    print "Number of target elements: " + str(len(atom_ids))

    n_pairs = (len(atom_ids)**2 - len(atom_ids) )/2
    print "Pairs of atoms: " + str(n_pairs)
    pairs_beg = 0
    pairs_end = n_pairs
    if opts.task_id > 0:
        n_pairs_in_task = (n_pairs + ( opts.n_div_job -1)) / opts.n_div_job
        pairs_beg = n_pairs_in_task * (opts.task_id-1)
        pairs_end = pairs_beg + n_pairs_in_task
    print "Pairs to be calculated in this job : " + str(pairs_beg) + " - " + str(pairs_end)
    task_atom_ids = check_target_atom_ids(atom_ids, pairs_beg, pairs_end)
    #print str(len(task_atom_ids))
    
    crd_files = None
    assign_files = None
    if opts.fn_crd_bin:
        crd_files, assign_files = check_files_bin(task_atom_ids, opts.fn_crd_bin,
                                                  opts.pref_assign, opts.suff_assign)
    else:
        crd_files, assign_files = check_files(task_atom_ids, opts.fn_gaussian,
                                              opts.fn_crd_bin,
                                              opts.pref_crd, opts.suff_crd, 
                                              opts.pref_assign, opts.suff_assign)

    print assign_files.keys()

    if opts.fn_o_mdcc and len(crd_files.keys())>0 and len(assign_files.keys()) > 0:
        print "mDCC"
        gaussian_info = read_gaussian_info(opts.fn_gaussian, task_atom_ids)
        cal_correlation(atom_ids, dim,
                        gaussian_info,
                        pairs_beg, pairs_end,  task_atom_ids,
                        opts.fn_o_mdcc,
                        opts.min_corr, opts.min_pi, 
                        opts.skip, crd_files, assign_files,
                        flg_crd_bin,
                        opts.flg_assign_bin,
                        opts.crd_tsv_columns, opts.crd_tsv_skip_header,
                        #opts.lbb_block_size, opts.lbb_b, opts.lbb_repeat,
                        time_range, opts.coef_mode)
    elif opts.fn_o_dcc and len(crd_files.keys())>0:
        print "DCC"
        cal_correlation_no_cluster(atom_ids, dim,
                                   pairs_beg, pairs_end, task_atom_ids,
                                   opts.fn_o_dcc,
                                   opts.min_corr,
                                   opts.skip,
                                   crd_files,
                                   flg_crd_bin,
                                   #opts.lbb_block_size,
                                   #opts.lbb_b, opts.lbb_repeat,
                                   opts.crd_tsv_columns, opts.crd_tsv_skip_header,
                                   time_range)
    elif opts.fn_o_rmsf and len(crd_files.keys())>0:
        print "RMSF"
        cal_rmsf(univ, atom_ids, 
                 opts.fn_o_rmsf,
                 opts.skip,
                 crd_files,
                 flg_crd_bin,
                 #opts.lbb_block_size,
                 #opts.lbb_b, opts.lbb_repeat,
                 opts.crd_tsv_columns, opts.crd_tsv_skip_header,
                 time_range, dim)
    else:
        print "--o-mdcc or --o-dcc or --o-rmsf must be specified"
    
    print "finished."
    #print "output_datatable"
    #output_corr_matrix(opts.fn_o_dat, univ, atom_ids, corr)


def check_target_atom_ids(atom_ids, pairs_beg, pairs_end):
    task_atom_ids = set()
    cnt = -1
    for a1 in atom_ids:
        for a2 in atom_ids:
            if a1 == a2: break
            cnt += 1
            if cnt < pairs_beg: continue
            if cnt >= pairs_end: continue
            task_atom_ids.add(a1)
            task_atom_ids.add(a2)
    print "check_target_atom_ids: " + str(cnt)
    print task_atom_ids
    return task_atom_ids

def select_atom_ids(str_select_id):
    ids = set()
    sp1 = str_select_id.split(",")
    for term in sp1:
        sp2 = term.split("-")
        if len(sp2) > 2:
            print "Error : the string for selection atom ids is wrong"
            print str_select_id
            sys.exit(1)
        for i in range(int(sp2[0]), int(sp2[1])+1):
            ids.add(i)
    return ids

def read_gaussian_info(fn_gaussian, atoms):
    gaussian_info = {}
    ## gaussian_info[atom_id] = (id, pi, mu, sigma)

    f = open(fn_gaussian)
    for line in f:
        terms = re.compile("\s+").split(line.strip())
        try:
            int(terms[0])
        except:
            continue
        gc_id = int(terms[0])
        atom = int(terms[1])
        pi = float(terms[2])
        mu = np.array([float(x) for x in terms[3:6]])
        sigma = np.array([[float(x) for x in terms[6:9]],
                             [float(x) for x in terms[9:12]],
                             [float(x) for x in terms[12:15]]])
        if len(atoms) == 0 or atom in atoms:
            if not atom in gaussian_info:
                gaussian_info[atom] = []
            gaussian_info[atom].append((gc_id, pi, mu, sigma))
    f.close()
    return gaussian_info

def check_files_bin(atom_ids, fn_crd_bin,
                    pref_assign, suff_assign):
    crd_files = {}
    assign_files = {}    
    if not os.path.isfile(fn_crd_bin):
        sys.stderr.write("Coordinate file " + fn_crd_bin + " is not found\n")
        #return None, None
    
    for aid in atom_ids:
        crd_files[aid] = fn_crd_bin
        if not pref_assign: continue
        fn_assign = pref_assign+str(aid)+suff_assign
        if not os.path.isfile(fn_assign):
            sys.stderr.write("Assign file " + fn_assign + " is not found\n")
            #return None, None
        assign_files[aid] = fn_assign

    return crd_files, assign_files
        
def check_files(atom_ids, fn_gaussian, fn_crd_bin, pref_crd, suff_crd,
                pref_assign, suff_assign):
    
    crd_files = {}
    assign_files = {}
    #print "The number of files " + len(os.listdir("."))
    for aid in atom_ids:
        crd_files[aid] = fn_crd_bin
        if not pref_crd: continue

        fn_crd = pref_crd+str(aid)+suff_crd
        if not os.path.isfile(fn_crd):
            sys.stderr.write("Coordinate file " + fn_crd + " is not found\n")
        crd_files[aid] = fn_crd

        fn_assign = pref_assign+str(aid)+suff_assign
        if not os.path.isfile(fn_assign):
            sys.stderr.write("Assign file " + fn_assign + " is not found\n")
            #return None, None
        assign_files[aid] = fn_assign

    flg_stop = False
    for atom_id in atom_ids:
        if not fn_crd_bin:
            if not atom_id in crd_files:
                sys.stderr.write("Coordinate file for atom " + str(atom_id) + " is not found\n")
                flg_stop = True
        else:
            crd_files[atom_id] = fn_crd_bin
        if pref_assign and not atom_id in assign_files:
            sys.stderr.write("Assignment file for atom " + str(atom_id) + " is not found\n")
            flg_stop = True
    #print crd_files
    #if flg_stop:
    #    sys.exit()

    return crd_files, assign_files
    
def get_atom_ids_from_select(univ, select):
    atom_ids = []
    for sel in select:
        atoms = univ.selectAtoms(sel)
        atom_ids.extend([atom.number for atom in atoms])
        #print sel + " : " + " ".join([str(atom.number) for atom in atoms])
        print sel + " : " + str(len([str(atom.number) for atom in atoms]))
    return sorted(set(atom_ids))

def cal_ave_crd(univ, atom_ids, time_range,
                skip):#, timestep, start_time):
    ## time_range is a tuple with 2 integers
    ##   that indicated first and last frames
    ## start_time means time of the first frame in
    ##   the trajectory file

    ## summation of all coordinates over tiem
    sum_crd = {}
    for atom_id in atom_ids:
        sum_crd[atom_id] = np.array([0.0,0.0,0.0])

    count = 0
    for n_frames, ts in enumerate(univ.trajectory):
        if n_frames % skip != 0: continue
        if n_frames % 1000 == 0: print n_frames
        if time_range[0] > n_frames: continue
        if time_range[1] > 0 and time_range[1] <= n_frames: break
        #time = ts.time
        #if timestep:
        #time = timestep * n_frames + start_time        
        #times.append(time)
        count += 1
        for atom_id in atom_ids:
            sum_crd[atom_id] += ts[atom_id]
            
    ave_crd = {}
    for atom_id, crd in sum_crd.items():
        ave_crd[atom_id] = crd/float(count)
    return ave_crd
            
def get_traj(univ, atom_id, skip, time_range):
    traj = [[],[],[]]
    for n_frames, ts in enumerate(univ.trajectory):
        if n_frames % skip != 0: continue
        if time_range[0] > n_frames: continue
        if time_range[1] > 0 and time_range[1] <= n_frames: break
        for i in [0,1,2]:
            traj[i].append(ts[atom_id][i])
    return np.array(traj[0]), np.array(traj[1]), np.array(traj[2])

def cal_fluct(univ, atom_id, skip,  time_range):
    traj = get_traj(univ, atom_id, skip,   time_range)
    ave = [np.mean(x) for x in traj]
    fluct = [ x - ave[i] for i, x in enumerate(traj) ]
    return fluct
    
def read_traj_crd(crd_files, atom_ids, columns, skip, skip_header, time_range=(0,-1)):
    traj_crd = {}
    for atom_id in atom_ids:
        traj = []
        f = open(crd_files[atom_id])
        print crd_files[atom_id]
        for i in range(skip_header): line_header = f.readline()
        
        traj = []
        for i in range(len(columns)): traj.append([])

        for time, line in enumerate(f):
            if time%skip!=0: continue
            terms = re.compile("\s+").split(line.strip())
            if time<time_range[0] or (time_range[1] > 0 and time>=time_range[1]): continue
            for i,x in enumerate(columns): traj[i].append(float(terms[x]))
            #traj.append(crd)
        traj_crd[atom_id] = np.array(traj)
        f.close()
        #print traj
    return traj_crd

def read_traj_pdf_bin(assign_file, read_from=0, read_to=-1):
    f = open(assign_file)
    flg_endian = ">"
    buf = f.read(4)
    magic = struct.unpack(">i", buf)[0]
    if magic != 1993:
        magic = struct.unpack("<i", buf)[0]
        flg_endian = "<"
        if magic != 1993:
            print "Invalid assign file format: the first four byte is not 1993: " + assign_file
            sys.exit()

    read_st = flg_endian+"i"
    n_col = struct.unpack(read_st, f.read(4))[0]
    n_frames = struct.unpack(read_st, f.read(4))[0]
    if not 0 <= read_to <= n_frames:
        read_to = n_frames        
    n_read_frames = read_to - read_from

    n_gauss = n_col/3
    pdf_array = []
    for i in range(n_gauss):
        f.seek(12+8*n_frames*(3*i+2) + 8*read_from)
        read_st = flg_endian+str(n_read_frames)+"d"
        buf = f.read(8*n_read_frames)        
        pdf = struct.unpack(read_st, buf)
        pdf_array.append(pdf)
    f.close()
    return np.array(pdf_array)

def read_traj_pdf(assign_file, time_range):
    ## assign_file
    ## Each line corresponds each time step
    ## The number of columns is #Gc * 3
    ## (repeat of [Maharanobis distance, pdf for Gc, pdf for GMM])
    ## 
    traj = []
    f = open(assign_file)
    n_read_lines = 0
    for i, line in enumerate(f):
        if time_range[0] > i: continue
        if time_range[1] >= 0 and i >=  time_range[1]: break
        #if n_read_lines == n_lines: break
        terms = re.compile("\s+").split(line.strip())
        tmp = np.array([float(terms[i+1]) for i in range(0, len(terms), 3)])
        traj.append(tmp)
        n_read_lines += 1
    f.close()
    return np.array(traj).T

def load_file_atom(atom_id, time_range, skip,
                   flg_crd_bin, flg_assign_bin,
                   crd_files, assign_files,
                   crd_tsv_columns, crd_tsv_skip_header,
                   trajreader=None):
    traj_crd = None
    if flg_crd_bin:
        traj_crd = trajreader.read_atom_crd(atom_id, read_from=time_range[0], read_to=time_range[1])
    else:
        traj_crd = read_traj_crd(crd_files, [atom_id], crd_tsv_columns, skip, crd_tsv_skip_header, time_range)[atom_id]

    traj_pdf = None
    if len(assign_files.keys()) > 0:
        if flg_assign_bin:
            traj_pdf = read_traj_pdf_bin(assign_files[atom_id], read_from=time_range[0], read_to=time_range[1])
        else: 
            #traj_pdf = (read_traj_pdf(assign_files[atom_id], len(traj_crd), time_range)[atom_id])
            traj_pdf = read_traj_pdf(assign_files[atom_id], time_range)
    
    return traj_crd, traj_pdf

def cal_correlation(atom_ids,  dim,
                    gaussian_info,
                    pairs_beg, pairs_end, task_atom_ids,
                    fn_out_gauss, min_corr, min_pi,
                    skip, crd_files, assign_files, 
                    flg_crd_bin, flg_assign_bin,
                    crd_tsv_columns, crd_tsv_skip_header,
                    #lbb_block_size, lbb_b, lbb_repeat,
                    time_range=(0,-1),
                    coef_mode="pkpl"):

    f_gauss = open(fn_out_gauss,"w")

    count = 0

    n_pairs = (len(atom_ids) * len(atom_ids) - len(atom_ids))/2
    cur_pair = -1
    trajreader = None
    if flg_crd_bin:
        print crd_files.values()[0]
        trajreader = kktrajtrans.TrajTReader(crd_files.values()[0])
        trajreader.open()
        trajreader.read_header()
    print atom_ids
    for i, atom_id1 in enumerate(atom_ids):
        if not atom_id1 in task_atom_ids:
            for atom_id2 in atom_ids:
                if atom_id1 == atom_id2: break
                cur_pair += 1
                if pairs_beg > cur_pair: continue
                if pairs_end <= cur_pair: break
            continue

        traj_crd1, traj_pdf1 = load_file_atom(atom_id1, time_range, skip,
                                              flg_crd_bin, flg_assign_bin,
                                              crd_files, assign_files,
                                              crd_tsv_columns, crd_tsv_skip_header,
                                              trajreader)
        
        valid_gc_in_atom1 = []
        for x, gc in enumerate(gaussian_info[atom_id1]):
            if gc[1] >= min_pi: valid_gc_in_atom1.append(x)
        for at in valid_gc_in_atom1:
            traj_pdf1[at] += EPS6

        n_frames = traj_crd1.shape[1]
        pdfsum1 = np.sum(traj_pdf1, axis=0)

        for atom_id2 in atom_ids:
            if atom_id1 == atom_id2: break
            cur_pair += 1
            if pairs_beg > cur_pair: continue
            if pairs_end <= cur_pair: break
            #if not atom_id1 in task_atom_ids: continue
            #if not atom_id2 in task_atom_ids: continue
            traj_crd2, traj_pdf2 = load_file_atom(atom_id2, time_range, skip,
                                                  flg_crd_bin, flg_assign_bin,
                                                  crd_files, assign_files,
                                                  crd_tsv_columns, crd_tsv_skip_header,
                                                  trajreader)

            valid_gc_in_atom2 = []
            for x, gc in enumerate(gaussian_info[atom_id2]):
                if gc[1] >= min_pi: valid_gc_in_atom2.append(x)
            traj_pdf2 = traj_pdf2[valid_gc_in_atom2] + EPS6

            pdfsum2 = np.sum(traj_pdf2, axis=0)
            
            def calc_terms():
                vgcindex = {}
                cross12 = []
                dot_xx1 = {}
                dot_xx2 = {}
                pi1 = {}
                pi2 = {}
                disp1 = {}
                for vgc1, gc_in_atom1 in enumerate(valid_gc_in_atom1):
                    g1 = gaussian_info[atom_id1][gc_in_atom1]
                    gc1 = g1[0]
                    mu1 = g1[2]
                    pi1[vgc1] = traj_pdf1[vgc1]/pdfsum1
                    disp1[vgc1] = np.array([traj_crd1[i] - mu1[i] for i in range(dim)]) + EPS6
                    dot_xx1[vgc1] = np.array([np.dot(x,x) for x in disp1[vgc1].T])

                for vgc2, gc_in_atom2 in enumerate(valid_gc_in_atom2):
                    g2 = gaussian_info[atom_id2][gc_in_atom2]
                    gc2 = g2[0]
                    mu2 = g2[2]
                    pi2[vgc2] = traj_pdf2[vgc2]/pdfsum2
                    disp2 = np.array([traj_crd2[i] - mu2[i] for i in range(dim)]) + EPS6
                    dot_xx2[vgc2] = np.array([np.dot(x,x) for x in disp2.T])
                    for vgc1, gc_in_atom1 in enumerate(valid_gc_in_atom1):
                        vgcindex[(vgc1, vgc2)] = len(cross12)
                        cross12.append(np.array([np.dot(d1,d2) for d1, d2 in zip(disp1[vgc1].T, disp2.T)]))
                return np.array(cross12), dot_xx1, dot_xx2, pi1, pi2, vgcindex
            cross12, dot_xx1, dot_xx2, pi1, pi2, vgcindex =  calc_terms()
            sum_cross12 = np.sum(np.abs(cross12), axis=0)

            for vgc1, gc_in_atom1 in enumerate(valid_gc_in_atom1):
                g1 = gaussian_info[atom_id1][gc_in_atom1]
                gc1 = g1[0]
                mu1 = g1[2]
                for vgc2, gc_in_atom2 in enumerate(valid_gc_in_atom2):
                    g2 = gaussian_info[atom_id2][gc_in_atom2]
                    gc2 = g2[0]
                    mu2 = g2[2]
                    pi12 = pi1[vgc1]*pi2[vgc2]
                    coef12 = pi12
                    coef1 = pi12
                    coef2 = pi12
                    weight_cross = 1.0
                    if coef_mode == "pkpk":
                        coef1 = pi1[vgc1]
                        coef2 = pi2[vgc2]
                    #elif coef_mode == "product_weight":
                    #    weight_cross = np.abs(cross12[vgcindex[(vgc1, vgc2)]]) / sum_cross12
                    #    coef12 = pi12 * weight_cross
                    #    coef1 = pi12 * weight_cross
                    #    coef2 = pi12 * weight_cross
                    cross = coef12 * cross12[vgcindex[(vgc1, vgc2)]]
                    cross_sum = np.sum(cross)
                    #dot1_w = pi12  * dot_xx1[vgc1]
                    dot1_w = coef1 * dot_xx1[vgc1]
                    dot1_w_sum = np.sum(dot1_w)
                    #dot2_w = pi12 * dot_xx2[vgc2]
                    dot2_w = coef2 * dot_xx2[vgc2]
                    dot2_w_sum = np.sum(dot2_w)
                    
                    corr = cross_sum / np.sqrt(dot1_w_sum * dot2_w_sum)

                    ### local block bootstrap
                    #    b_corr_sum = 0.0
                    #    b_cross_sum = 0.0
                    #    b_dot1_w_sum = 0.0
                    #    b_dot2_w_sum = 0.0
                    #    b_corr_sum2 = 0.0
                    #    b_cross_sum2 = 0.0
                    #    b_dot1_w_sum2 = 0.0
                    #    b_dot2_w_sum2 = 0.0
                    #    
                    #    for i in xrange(lbb_repeat):
                    #        series = lbb.gen_series(n_frames, lbb_block_size, lbb_b)
                    #        cur_dot1_w = np.sum(dot1_w[series])
                    #        cur_dot2_w = np.sum(dot2_w[series])
                    #        cur_cross  = np.sum(cross[series])
                    #        cur_corr   = cur_cross / np.sqrt(cur_dot1_w*cur_dot2_w)
                    #        b_dot1_w_sum += cur_dot1_w
                    #        b_dot2_w_sum += cur_dot2_w
                    #        b_cross_sum  += cur_cross
                    #        b_corr_sum += cur_corr
                    #        b_dot1_w_sum2 += cur_dot1_w*cur_dot1_w
                    #        b_dot2_w_sum2 += cur_dot2_w*cur_dot2_w
                    #        b_cross_sum2  += cur_cross*cur_cross
                    #        b_corr_sum2   += cur_corr*cur_corr

                    #    mean_corr = b_corr_sum/float(lbb_repeat)
                    #    sd_corr = np.sqrt(b_corr_sum2/float(lbb_repeat) - mean_corr*mean_corr)
                    #    mean_cross = b_cross_sum/float(lbb_repeat)
                    #    sd_cross = np.sqrt(b_cross_sum2/float(lbb_repeat) - mean_cross*mean_cross)
                        
                    d_mu = (mu1 - mu2)
                    dist_mu = np.sqrt(np.sum(d_mu*d_mu))
                    
                    #Pi_coef = np.sum(pi)/pi12sum
                    coef_corr = np.sum(coef12) / (np.sqrt(np.sum(coef1)*np.sum(coef2)))

                    if np.fabs(corr) >= min_corr:
                        line = "\t".join([str(x) for x in [gc2, gc1, atom_id2, atom_id1, corr,
                                                           np.average(pi12), dist_mu,
                                                           ##coef_corr, 
                                                           #np.average(cross),
                                                           ##np.average(dot1_w), np.average(dot2_w)
                                                           #np.average(coef12), np.average(coef1),  np.average(coef2)
                                                           ] ] )
                                                           #, mean_corr, sd_corr, mean_cross, sd_cross]])
                                                           #pi_coef, pi_coef * corr] ] )
                        f_gauss.write(line + "\n")

    f_gauss.close()
    if flg_crd_bin: trajreader.close()        
    return 

def cal_correlation_no_cluster(atom_ids, dim,
                               pairs_beg, pairs_end, task_atom_ids,
                               fn_out_dcc, min_corr,
                               skip, crd_files, flg_crd_bin,
                               #lbb_block_size, lbb_b, lbb_repeat,
                               crd_tsv_columns, crd_tsv_skip_header,
                               time_range=(0,-1)):

    f_dcc = open(fn_out_dcc,"w")

    count = 0

    dot_xx = {}

    n_pairs = (len(atom_ids) * len(atom_ids) - len(atom_ids))/2
    cur_pair = -1
    trajreader = None
    if flg_crd_bin:
        print crd_files.values()[0]
        trajreader = kktrajtrans.TrajTReader(crd_files.values()[0])
        trajreader.open()
        trajreader.read_header()
    dummy = {}
    print atom_ids
    for i, atom_id1 in enumerate(atom_ids):
        #print "atom : " + str(atom_id1)
        #if not atom_id1 in task_atom_ids: continue
        #print " calc DCC " + str(i) + " " + str(cur_pair) + " " + str(pairs_beg) +"-"+str(pairs_end)
        traj_crd1, traj_pdf1 = load_file_atom(atom_id1, time_range, skip,
                                              flg_crd_bin, False,
                                              crd_files, dummy,
                                              crd_tsv_columns, crd_tsv_skip_header,
                                              trajreader)

        n_frames = traj_crd1.shape[1]
        mu1 = None
        
        #if center_mode=="origin":
        #    mu1 = np.array([0.0, 0.0, 0.0])
        #else:
        mu1 = np.average(traj_crd1, axis=1)
        disp1 = np.array([traj_crd1[i] - mu1[i] for i in range(dim)])
        dot_xx1 = np.array([np.dot(x,x) for x in disp1.T])
        for atom_id2 in atom_ids:
            if atom_id1 == atom_id2: break
            cur_pair += 1
            if pairs_beg > cur_pair: continue
            if pairs_end <= cur_pair: continue
            if not atom_id1 in task_atom_ids: continue
            if not atom_id2 in task_atom_ids: continue
            traj_crd2, traj_pdf2 = load_file_atom(atom_id2, time_range, skip,
                                                  flg_crd_bin, False,
                                                  crd_files, dummy,
                                                  crd_tsv_columns, crd_tsv_skip_header,
                                                  trajreader)
            #print str(atom_id1) + " " + str(atom_id2)
            mu2 = None
            #if center_mode=="origin":
            #    mu2 = np.array([0.0, 0.0, 0.0])
            #else:
            mu2 = np.average(traj_crd2, axis=1)
            disp2 = np.array([traj_crd2[i] - mu2[i] for i in range(dim)])
            dot_xx2 = np.array([np.dot(x,x) for x in disp2.T])

            cross = np.array([np.dot(d1,d2) for d1, d2 in zip(disp1.T, disp2.T)])
            cross_sum = np.sum(cross)
                    #cross = np.sum([pi[frame] * np.dot(d1, disp2[frame]) for frame, d1 in enumerate(disp1.T)])
            dot1_w = dot_xx1
            dot1_w_sum = np.sum(dot1_w)
            dot2_w = dot_xx2
            dot2_w_sum = np.sum(dot2_w)
                    
            corr = cross_sum / np.sqrt(dot1_w_sum * dot2_w_sum)
            mean_corr = 0.0
            mean_cross = 0.0
            sd_corr = 0.0
            sd_cross = 0.0

            d_mu = (mu1 - mu2)
            dist_mu = np.sqrt(np.sum(d_mu*d_mu))

            if np.fabs(corr) >= min_corr or np.fabs(mean_corr) >= min_corr:
                line = "\t".join([str(x) for x in [atom_id2, atom_id1, corr, dist_mu]])
                f_dcc.write(line + "\n")
                
    f_dcc.close()
    if flg_crd_bin: trajreader.close()        
    return 

def cal_rmsf(univ, atom_ids, 
             fn_o_rmsf,
             skip,
             crd_files,
             flg_crd_bin,
             #lbb_block_size,
             #lbb_b, lbb_repeat,
             crd_tsv_columns, crd_tsv_skip_header,
             time_range,
             dim):

    f_rmsf = open(fn_o_rmsf,"w")
    n_pairs = (len(atom_ids) * len(atom_ids) - len(atom_ids))/2
    cur_pair = -1
    trajreader = None

    if flg_crd_bin:
        print crd_files.values()[0]
        trajreader = kktrajtrans.TrajTReader(crd_files.values()[0])
        trajreader.open()
        trajreader.read_header()
    dummy = {}
    for atom_id in sorted(atom_ids):
        #print atom_id
        traj_crd, traj_pdf = load_file_atom(atom_id, time_range, skip,
                                            flg_crd_bin, False,
                                            crd_files, dummy,
                                            crd_tsv_columns, crd_tsv_skip_header,
                                            trajreader)
        n_frames = traj_crd.shape[1]

        mu = np.average(traj_crd, axis=1)
        disp = np.array([traj_crd[i] - mu[i] for i in range(dim)])
        #disp2 = np.array([x*x for x in disp.T])
        #disp2 = np.array( np.sum([x[i]*x[i] for i in [0,1,2]])  )
        disp2 = disp * disp
        mean_rmsf = 0.0
        sd_rmsf = 0.0
        #rmsf = np.sqrt(np.average(disp2))
        rmsf = np.sqrt(np.sum(disp2)/float(n_frames))
        
        #if lbb_block_size and lbb_b:
        #    ## local block bootstrap
        #    b_rmsf_sum = 0.0
        #    b_rmsf_sum2 = 0.0
        #    for i in xrange(lbb_repeat):
        #        series = lbb.gen_series(n_frames, lbb_block_size, lbb_b)
        #        cur_msf = np.sum(disp2[series])/float(n_frames)
        #        b_rmsf_sum += np.sqrt(cur_msf)
        #        b_rmsf_sum2 += cur_msf
        
        #    mean_rmsf = b_rmsf_sum/float(lbb_repeat)
        #    sd_rmsf = np.sqrt(b_rmsf_sum2/float(lbb_repeat) - mean_rmsf*mean_rmsf)

        line = "\t".join([str(x) for x in [atom_id, univ.atoms[atom_id].name,
                                           univ.atoms[atom_id].resnum,
                                           univ.atoms[atom_id].resname,
                                           rmsf]]) #, mean_rmsf, sd_rmsf]])
        f_rmsf.write(line + "\n")

    f_rmsf.close()
    return 
            
def output_corr_matrix(fn_o_dat, univ, atom_ids, corr):
    header = "\t".join([str(x)+"."+univ.atoms[x].name.replace("'",".") for x in atom_ids])
    f=open(fn_o_dat,"w")
    f.writer(header+"\n")
    for atom_id1 in atom_ids:
        line = str(atom_id1)+"."+univ.atoms[atom_id1].name.replace("'",".")
        for atom_id2 in atom_ids:
            pair = (atom_id2, atom_id1)
            if atom_id2 > atom_id1: pair = (atom_id1, atom_id2)
            if pair[0] == pair[1]: line+="\t1.0"
            else:
                line += "\t"+str(corr[pair])
        f.write(line + "\n")
    f.close()
    return 

if __name__ == "__main__":
    _main()
