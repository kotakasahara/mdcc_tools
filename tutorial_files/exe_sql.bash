#!/bin/bash
#$ -S /bin/bash
#$ -cwd

dbfile="data.db"
sqlite3 ${dbfile} < sqlquery_01.sql
sqlite3 ${dbfile} < sqlquery_02.sql
sqlite3 ${dbfile} < sqlquery_03.sql

dbfile_nw="data_nw.db"
sqlite3 ${dbfile_nw} < sqlquery_11.sql
sqlite3 ${dbfile_nw} < sqlquery_02.sql
sqlite3 ${dbfile_nw} < sqlquery_13.sql

echo "res_num1.int res_num2.int atom_id1.int atom_id2.int gauss_id1.int gauss_id2.int correlation.float coef.float dist.float corr_dcc.float dist_ave.float dist_sd.float dist_min.float dist_max.float dist_init.float" > tmp.txt
cat tmp.txt res_res.txt > res_res_h.txt
cat tmp.txt res_res_d5_c50.txt > res_res_d5_c50_h.txt
mv res_res_h.txt res_res.txt
mv res_res_d5_c50_h.txt res_res_d5_c50.txt

echo "atom_id1.int atom_id2.int gauss_id1.int gauss_id2.int correlation.float coef.float dist.float corr_dcc.float dist_ave.float dist_sd.float dist_min.float dist_max.float dist_init.float" > tmp.txt
cat tmp.txt atom_atom.txt > atom_atom_h.txt
cat tmp.txt atom_atom_d5_c50.txt > atom_atom_d5_c50_h.txt
mv atom_atom_h.txt atom_atom.txt
mv atom_atom_d5_c50_h.txt atom_atom_d5_c50.txt

python2.7 ${MDCCTOOLS}/bin/nx_centrality.py \
  -i res_res_d5_c50.txt \
  --i-elem ../ref_res.txt \
  --key-elem res_id \
  -o res_cent_btw_d5_c50.txt \
  --btw 
