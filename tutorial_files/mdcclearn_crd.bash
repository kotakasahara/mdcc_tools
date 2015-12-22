#!/bin/bash
#$ -S /bin/bash
#$ -cwd

n_cal=${1}
id_cal=${2}

fn_mdcclearn_conf="mdcclearn_conf_template.txt"
fn_mdcclearn_sh="mdcclearn_${id_cal}.bash"
fn_pdb="reference.pdb"

echo "${id_cal} / ${n_cal}"

python ${MDCCTOOLS}/bin/run_mdcc_tool.py \
    --mode mdcclearn \
    --pdb ${fn_pdb} \
    --select "not type H" \
    --n-div ${n_cal} \
    --task-id ${id_cal} \
    --mdcc-conf ${fn_mdcclearn_conf} \
    --mdcc-bin "${MDCCTOOLS}/bin/mdcc_learn" \
    --fn-sh mdcclearn_bash/${fn_mdcclearn_sh} 
    
cd mdcclearn_bash
bash ./${fn_mdcclearn_sh} > log_${id_cal}.txt
