#!/bin/bash
#$ -S /bin/bash
#$ -cwd

n_cal=${1}
id_cal=${2}

fn_mdccassign_conf="mdccassign_conf_template.txt"
fn_mdccassign_sh="mdccassign_${id_cal}.bash"
fn_pdb="reference.pdb"

echo "${id_cal} / ${n_cal}"

python2.7 ${MDCCTOOLS}/bin/run_mdcc_tool.py \
    --mode mdccassign \
    --pdb ${fn_pdb} \
    --select "not type H" \
    --n-div ${n_cal} \
    --task-id ${id_cal} \
    --mdcc-conf ${fn_mdccassign_conf} \
    --mdcc-bin "${MDCCTOOLS}/bin/mdcc_assign" \
    --fn-sh mdccassign_bash/${fn_mdccassign_sh} 
    
cd mdccassign_bash
bash ./${fn_mdccassign_sh}
