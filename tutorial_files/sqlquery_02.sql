.mode tabs

CREATE TABLE corr_dcc(
  atom_id1 int unsigned,   atom_id2 int unsigned,
  correlation float,    dist float);

.import ../corr_dcc.txt  corr_dcc

CREATE INDEX corrnc_atoms ON corr_dcc (atom_id1, atom_id2);
CREATE INDEX corrnc_dist ON corr_dcc (dist);
CREATE INDEX corrnc_corr ON corr_dcc (correlation);

CREATE TABLE corr_atom_nc AS SELECT c.*, nc.correlation AS corr_dcc FROM corr_atom AS c INNER JOIN corr_dcc AS nc ON c.atom_id1=nc.atom_id1 AND c.atom_id2=nc.atom_id2;

/* res-res */

CREATE INDEX corrnc_atom1 ON corr_dcc (atom_id1);
CREATE TEMPORARY TABLE tmp_corr_res1 AS SELECT corr_dcc.*, atoms.res_num AS res_num1 FROM corr_dcc INNER JOIN atoms ON corr_dcc.atom_id1=atoms.atom_id;
CREATE INDEX tmp_corrnc_res1 ON tmp_corr_res1 (atom_id2);
CREATE TEMPORARY TABLE tmp_corr_res2 AS SELECT tmp_corr_res1.*, atoms.res_num AS res_num2 FROM tmp_corr_res1 INNER JOIN atoms ON tmp_corr_res1.atom_id2=atoms.atom_id;

CREATE INDEX tmp_corrnc_dist ON tmp_corr_res2 (correlation);
CREATE TEMPORARY TABLE tmp_maxcorr_res AS SELECT res_num1, res_num2, max(correlation) AS correlation FROM tmp_corr_res2 GROUP BY res_num1, res_num2;
CREATE INDEX tmp_maxcorr_rrcor ON tmp_maxcorr_res (res_num1, res_num2, correlation);

CREATE TEMPORARY TABLE tmp_corr_res3 AS SELECT
  tmp2.res_num1, tmp2.res_num2,
  tmp2.atom_id1, tmp2.atom_id2,
  tmp2.correlation,  tmp2.dist
  FROM tmp_corr_res2 AS tmp2 INNER JOIN tmp_maxcorr_res AS tmpmax ON tmp2.res_num1 = tmpmax.res_num1 AND tmp2.res_num2 = tmpmax.res_num2 AND tmp2.correlation = tmpmax.correlation; 

CREATE TABLE corr_dcc_res AS SELECT DISTINCT * FROM tmp_corr_res3;
CREATE INDEX corr_dcc_res_id ON corr_dcc_res (res_num1, res_num2);
CREATE INDEX corr_res_id ON corr_res (res_num1, res_num2);

CREATE TABLE corr_res_nc AS SELECT c.*, nc.correlation AS corr_dcc FROM corr_res AS c INNER JOIN corr_dcc_res AS nc ON c.res_num1=nc.res_num1 AND c.res_num2=nc.res_num2;

