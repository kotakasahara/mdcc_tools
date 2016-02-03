#!/bin/bash
#$ -S /bin/bash
#$ -cwd

n_cal=${1}
id_cal=${2}

python2.7 ${MDCCTOOLS}/bin/cal_mdcc.py \
    --fn-ref reference.pdb \
    --gaussian crd_mdcclearn_gauss.txt \
    --pref-assign mdccassign_bash/assign/assign.dat. \
    --suff-assign "" \
    --fn-crd-bin  traj.trrmdcc \
    --select "not type H" \
    --o-mdcc mdcc/corr.txt.${id_cal} \
    --skip 1 \
    --n-div ${n_cal} \
    --min-corr 0.0 \
    --task-id ${id_cal} \
    --assign-binary \
    --frame-begin 10
