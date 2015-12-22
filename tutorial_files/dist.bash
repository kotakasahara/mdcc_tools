#!/bin/bash
#$ -S /bin/bash
#$ -cwd

source ~/.zshrc

python2.7  ${MDCCTOOLS}/bin/interatomic_distances.py \
  -i ../traj.mdcctrr \
  --init ../reference.pdb \
  --start-time 10000.0 \
  --atom-select "not type H and not name NA and not name CL" \
  --frame-begin 100 --frame-end -1\
  --n-div $1  --task-id ${SGE_TASK_ID} \
  -o dist.txt.${SGE_TASK_ID}
