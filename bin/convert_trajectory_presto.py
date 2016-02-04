#!/usr/bin/python2.6

from optparse import OptionParser
import numpy
import sys
import kkpdb
import kkpresto_crd as kkcrd
import kkgro_trr
import kkmdconf
import re
import copy

#time_factor = 0.1

def get_options():
    p = OptionParser()
    p.add_option('--i-conf', dest='fn_mdconf',
                 help="file name for mcdonf")
    p.add_option('--i-pdb', dest='fn_pdb',
                 help="file name for a presto initial structure")
    p.add_option('--remove-res', dest='rm_res',
                 action="append",
                 help="resname to be removed")
    p.add_option('--i-crd', dest='fn_crd',
                 help="file name for a presto coordinates")
    p.add_option('--i-vel', dest='fn_vel',
                 help="file name for a presto velocities")
    p.add_option('-o', dest='fn_trr',
                 help="file name for the output gromacs .trr file")
    p.add_option('--lx', dest='cell_x',
                 type = "float",
                 help="fixed cell length of x axis")
    p.add_option('--ly', dest='cell_y',
                 type = "float",
                 help="fixed cell length of y axis")
    p.add_option('--lz', dest='cell_z',
                 type = "float",
                 help="fixed cell length of z axis")
    p.add_option('-e', dest='end_frame',
                 type = "int",
                 help="last frame")

    p.add_option('--ow-frame-begin', dest='ow_frame_begin',
                 type = "int", default=-1,
                 help="The first frame, for over writing")
    p.add_option('--ow-time-step', dest='ow_time_step',
                 type = "float", default=-1.0,
                 help="time step")
    p.add_option('--ow-frame-step', dest='ow_frame_step',
                 type = "int", default=-1,
                 help="time step")

    p.add_option('--gro-presto', dest='flg_gro_presto',
                 action="store_true",
                 help="Flag for converting gro .trr file to presto format")
    
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args


def _main():
    opts, args = get_options()
    print "trajconv_presto_gro"
    ## set box information
    box = None
    if opts.cell_x and opts.cell_y and opts.cell_z:
        box = numpy.array([opts.cell_x, opts.cell_y, opts.cell_z])
    fn_pdb = opts.fn_pdb

    
    ## if mdconf file is input, read it
    mdconf = None
    if opts.fn_mdconf:
        mdconf = kkmdconf.MDConfReader(opts.fn_mdconf).read_conf()
        if not fn_pdb: fn_pdb = mdconf.psy_init_structure

    ## a list of atom ids which will be removed
    rm_atom_ids = set()
    if fn_pdb and opts.rm_res:
        struct = None
        struct = kkpdb.PDBReader(fn_pdb).read_model()
        rm_atom_ids = set(struct.get_atom_index_by_res_name(opts.rm_res))
        print "Ignoring atoms: " + str(len(rm_atom_ids))

    ## presto => gro
    if not opts.flg_gro_presto:
        if opts.fn_crd and opts.fn_trr:
            convert(opts.fn_trr, opts.fn_crd, opts.fn_vel, box=box, rm_atom_ids=rm_atom_ids,
                    frame_begin=opts.ow_frame_begin,
                    time_step=opts.ow_time_step,
                    frame_step=opts.ow_frame_step)  
        elif mdconf:
            cur_frame = opts.ow_frame_begin
            #print mdconf.psy_traj_crd
            for i, fn_crd in enumerate(mdconf.psy_traj_crd):
                fn_vel = None
                if len(mdconf.psy_traj_vel) != 0:
                    fn_vel = mdconf.psy_traj_vel[i]
                m = re.compile("(\d+)\.").search(fn_crd)
                run_id_num = i
                if m:
                    run_id_num = int(m.group(1))
                st = "%0"+str(mdconf.n_digit)+"d"
                run_id = st%run_id_num
                fn_trr = opts.fn_trr + run_id + ".trr"
                if fn_vel:
                    print fn_crd + ", " + fn_vel + " : " + fn_trr
                else:
                    print fn_crd + " : " + fn_trr
                read_frame = convert(fn_trr, fn_crd, fn_vel, box=box, rm_atom_ids=rm_atom_ids,
                                     frame_begin=cur_frame,
                                     time_step=opts.ow_time_step,
                                     frame_step=opts.ow_frame_step)
                cur_frame += read_frame * optw.ow_frame_step
    else:
        convert_gp(opts.fn_trr, opts.fn_crd,
                   rm_atom_ids=rm_atom_ids,
                   frame_begin=opts.ow_frame_begin,
                   time_step=opts.ow_time_step,
                   frame_step=opts.ow_frame_step)  
        
            
def convert(fn_trr, fn_crd, fn_vel=None, box=None, rm_atom_ids=set(),
            frame_begin=-1, time_step=-1.0, frame_step=-1):    
    cur_frame = frame_begin
    ## box : numpy.array([float, float, float])
    crd_reader = kkcrd.PrestoCrdReader(fn_crd)
    vel_reader = kkcrd.PrestoCrdReader(fn_vel)
    trr_writer = kkgro_trr.GroTrrWriter(fn_trr)

    crd_reader.open()
    if fn_vel: vel_reader.open()
    trr_writer.open()

    frame_crd = crd_reader.read_next_frame()
    read_frame = 1
    while frame_crd:
        #print frame_crd.step
        frame_vel_crds = None
        if fn_vel:
            frame_vel_crds = vel_reader.read_next_frame().crds
            #frame_vel_crds = scale_crds(frame_vel_crds)
            if not frame_vel:  break
            
        #frame_crd_crds = scale_crds(frame_crd.crds)
        w_time = frame_crd.time
        w_frame = frame_crd.step
        if time_step > 0:
            w_frame = cur_frame
            w_time = cur_frame * time_step
        #print "time: " + str(w_time) + " frame:" + str(w_frame)
        trr_writer.write_frame(box,
                               frame_crd.crds,
                               frame_vel_crds,
                               None,
                               w_frame, w_time,
                               bin = False,
                               prec = kkgro_trr.GroTrrWriter.PREC_SINGLE,
                               ignore = rm_atom_ids,
                               scale_length = 0.1)
        frame_crd = crd_reader.read_next_frame()
        cur_frame += frame_step
        read_frame += 1
    trr_writer.close()
    if fn_vel: vel_reader.close()    
    crd_reader.close()    

    return read_frame

def convert_gp(fn_trr, fn_crd, rm_atom_ids=set(),
               frame_begin=-1, time_step=-1.0, frame_step=-1):
    cur_frame = frame_begin
    ## box : numpy.array([float, float, float])
    trr_reader = kkgro_trr.GroTrrReader(fn_crd)
    crd_writer = kkcrd.PrestoCrdWriter(fn_trr)
    
    trr_reader.open()
    crd_writer.open()

    frame_crd = trr_reader.read_next_frame()
    read_frame = 1
    w_time = frame_crd.time
    w_frame = frame_crd.step
    print "time: " + str(w_time) + " frame:" + str(w_frame)
    while frame_crd:
        #print frame_crd.step
        #frame_crd_crds = scale_crds(frame_crd.crds)
        w_time = frame_crd.time
        w_frame = frame_crd.step
        if time_step > 0:
            w_frame = cur_frame
            w_time = cur_frame * time_step
        
        frame_crd.crds *= 10.0
        crd_writer.write_frame(frame_crd.crds,
                               w_frame, w_time,
                               bin = False,
                               ignore = rm_atom_ids)
        frame_crd = trr_reader.read_next_frame()
        cur_frame += frame_step
        read_frame += 1
    print "time: " + str(w_time) + " frame:" + str(w_frame)
    crd_writer.close()
    trr_reader.close()    
    return read_frame

def scale_crds(crds, scale = 0.1):
    for crd in crds: crd *= scale
    return crds

if __name__ == "__main__":
    _main()
