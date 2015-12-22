.mode tabs

CREATE TABLE dist(
  atom_id1 int unsigned,   atom_id2 int unsigned,
  dist_ave float, dist_sd float, dist_min float, dist_max float);

.import ../dist.txt dist

CREATE INDEX dist_atoms ON dist (atom_id1, atom_id2);

CREATE TABLE dist_init(
  atom_id1 int unsigned,   atom_id2 int unsigned,
  dist_init float);

.import ../dist_init.txt dist_init

CREATE INDEX dist_init_atoms ON dist_init (atom_id1, atom_id2);

CREATE INDEX corr_atom_nc_id ON corr_atom_nc (atom_id1, atom_id2);

CREATE TEMPORARY TABLE atom_atom_pre AS SELECT 
  cr.atom_id1,  cr.atom_id2, cr.gauss_id1,  cr.gauss_id2, 
  cr.correlation,  cr.coef,  cr.dist,  
  cr.corr_dcc,
  dist.dist_ave, dist.dist_sd, dist.dist_min, dist.dist_max 
  FROM corr_atom_nc AS cr INNER JOIN dist ON cr.atom_id1=dist.atom_id1 AND cr.atom_id2=dist.atom_id2;

CREATE INDEX tmp_aa ON atom_atom_pre (atom_id1, atom_id2);

CREATE TABLE atom_atom AS SELECT 
  cr.atom_id1,  cr.atom_id2, cr.gauss_id1,  cr.gauss_id2, 
  cr.correlation,  cr.coef,  cr.dist, 
  cr.corr_dcc,
  cr.dist_ave, cr.dist_sd, cr.dist_min, cr.dist_max,
  dist_init.dist_init
  FROM atom_atom_pre AS cr INNER JOIN dist_init ON cr.atom_id1=dist_init.atom_id1 AND cr.atom_id2=dist_init.atom_id2;

/* res-res */
.mode tabs
CREATE TEMPORARY TABLE res_res_pre AS SELECT
  cr.*, dist.dist_ave, dist.dist_sd, dist.dist_min, dist.dist_max
  FROM corr_res_nc AS cr 
  INNER JOIN dist ON cr.atom_id1=dist.atom_id1 AND cr.atom_id2=dist.atom_id2;
CREATE TABLE res_res AS SELECT
  cr.*, dist_init.dist_init
  FROM res_res_pre AS cr 
  INNER JOIN dist_init ON cr.atom_id1=dist_init.atom_id1 AND cr.atom_id2=dist_init.atom_id2;

.output res_res.txt
SELECT * FROM res_res;

.output atom_atom.txt
SELECT * FROM atom_atom;

