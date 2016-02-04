#!/bin/bash
#$ -S /bin/bash
#$ -cwd

n_cal=${1}
id_cal=${2}

python2.7 ${MDCCTOOLS}/bin/cal_mdcc.py \
    --fn-ref reference_cano.pdb \
    --fn-crd-bin  traj.trrmdcc \
    --select "not type H" \
    --o-dcc dcc/corr.txt.${id_cal} \
    --skip 1 \
    --n-div ${n_cal} \
    --min-corr 0.0 \
    --task-id ${id_cal} \
    --frame-begin 10
