#!/bin/bash
#$ -S /bin/bash
#$ -cwd

python2.7 ${MDCCTOOLS}/bin/interatomic_distances_pdb.py  \
  --i-pdb reference.pdb \
  --select "not type H" \
  --o-dist dist_init.txt
