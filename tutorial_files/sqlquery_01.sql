.mode tabs

CREATE TABLE corr(
  gauss_id1 int unsigned,  gauss_id2 int unsigned,
  atom_id1 int unsigned,   atom_id2 int unsigned,
  correlation float, coef float,
  dist float);

.import ../corr_mdcc.txt  corr

CREATE INDEX corr_atoms ON corr (atom_id1, atom_id2);
/* CREATE INDEX corr_dist ON corr (dist); */
/* CREATE INDEX corr_corr ON corr (correlation); */

CREATE TABLE atoms(
  atom_id int unsigned,    atom_name  char(32),
  res_name char(32),       res_num int unsigned,
  res_id int unsigned,     res_id_auth int unsigned,
  chain_id char(32),       segment char(32) , seg_type char(32)  );

.import ../ref_atoms_woh.txt atoms
CREATE INDEX atom_id ON atoms(atom_id); 

CREATE TABLE res(
  res_num int unsigned,     res_id int unsigned,
  res_id_auth int unsigned, res_name char(32),
  first_atom_id int unsigned, last_atom_id int unsigned,
  res_label char(32), chain_id char(32),
  segment char(32), seg_type char(32)
);
.import ../ref_res_woh.txt res


CREATE TABLE crd_gauss(
  gauss_id int unsigned,   atom_id int unsigned,
  pi float unsigned,    mu1 float,     mu2 float,     mu3 float,
  sig11 float,     sig12 float,   sig13 float,
  sig21 float,     sig22 float,   sig23 float,
  sig31 float,     sig32 float,   sig33 float);

.import ../crd_mdcclearn_gauss.txt crd_gauss

CREATE INDEX crd_gauss_id ON crd_gauss (gauss_id);
CREATE INDEX crd_gauss_atomid ON crd_gauss (atom_id);

CREATE TABLE gauss_atominfo AS SELECT * FROM crd_gauss 
INNER JOIN atoms ON crd_gauss.atom_id = atoms.atom_id;

CREATE INDEX crdga_atomid ON gauss_atominfo(atom_id);
CREATE INDEX crdga_gaussid ON gauss_atominfo(gauss_id);

/* atom-atom corr */

/* find max coef g-g pair with coef>=0.1 in a-a pair */
CREATE TABLE corr_atom AS SELECT 
  gauss_id1, gauss_id2, atom_id1, atom_id2,
  max(correlation) AS correlation, coef, dist 
  FROM corr WHERE coef >= 0.1  GROUP BY atom_id1, atom_id2;

/* res-res corr */

CREATE INDEX corr_atom_id12 ON corr_atom (atom_id1, atom_id2);
CREATE INDEX corr_atom_id1 ON corr_atom (atom_id1);
CREATE TEMPORARY TABLE tmp_corr_res1
  AS SELECT corr_atom.*, atoms.res_num AS res_num1, atoms.res_id_auth AS res_id1, atoms.res_name AS res_name1
  FROM corr_atom INNER JOIN atoms ON corr_atom.atom_id1=atoms.atom_id;

CREATE INDEX tmp_corr_res1_atom2 ON tmp_corr_res1 (atom_id2);

CREATE TEMPORARY TABLE tmp_corr_res2
  AS SELECT tmp_corr_res1.*, atoms.res_num AS res_num2, atoms.res_id_auth AS res_id2, atoms.res_name AS res_name2
  FROM tmp_corr_res1 INNER JOIN atoms ON tmp_corr_res1.atom_id2=atoms.atom_id;

CREATE INDEX tmp_corr_res2_corr ON tmp_corr_res2 (correlation);

CREATE TEMPORARY TABLE tmp_maxcorr_res AS SELECT res_num1, res_num2, max(correlation) AS correlation 
  FROM tmp_corr_res2 GROUP BY res_num1, res_num2;

CREATE INDEX tmp_corr_res2_rr_corr ON tmp_corr_res2 (res_num1, res_num2, correlation);
CREATE INDEX tmp_maxcorr_rr_corr ON tmp_maxcorr_res (res_num1, res_num2, correlation);

CREATE TEMPORARY TABLE tmp_corr_res3 AS SELECT
  tmp2.res_num1, tmp2.res_num2,
  tmp2.atom_id1, tmp2.atom_id2,
  tmp2.gauss_id1, tmp2.gauss_id2,
  tmp2.correlation, tmp2.coef, tmp2.dist
 FROM tmp_corr_res2 AS tmp2 INNER JOIN tmp_maxcorr_res AS tmpmax 
  ON tmp2.res_num1 = tmpmax.res_num1 AND 
     tmp2.res_num2 = tmpmax.res_num2 AND
     tmp2.correlation = tmpmax.correlation;

CREATE TABLE corr_res AS SELECT DISTINCT * FROM tmp_corr_res3;

