#!/bin/bash
#$ -S /bin/bash
#$ -cwd

n_cal=${1}
id_cal=${2}


python2.7  ${MDCCTOOLS}/bin/interatomic_distances.py \
  -i traj.trrmdcc \
  --init reference.pdb \
  --atom-select "not type H" \
  --n-div ${n_cal}  --task-id ${id_cal} \
  --frame-begin 10 \
  -o dist/dist.txt.${id_cal}
