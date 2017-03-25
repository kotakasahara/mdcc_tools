#include "PliGauss.h"

PliGauss::PliGauss(){}

PliGauss::~PliGauss(){

}

int PliGauss::setup(int argn, char* argv[]){
  string fn_cfg;
  if(argn<2){
    cerr << "Usage: has [mode]"<<endl;
    cerr << "------------------------------------"<<endl;
    exit(1);
  }
  cout << "conf.setall\n";
  cfg.set_all(argn,argv);
  if(cfg.fn_cfg!=string())
    cfg.set_all(Read(cfg.fn_cfg).load_config());
  cout <<"/setup\n";
  

  return 0;
}

int PliGauss::mainStream(){
  cout<<"mainstream\n";
  switch(cfg.mode){
  case M_TEST: test(); break;
  case M_ASSIGN_TTI: assign_tti(); break;
  case M_ASSIGN_TAI: assign_tai(); break;
  case M_ASSIGN_TABLE: assign_datatable(); break;
  case M_ASSIGN_TRAJTRANS: assign_trajtrans(); break;
  default:
    cout <<"mode ga shitei sareteimasen.\n";
    break;
  }
  cout <<"finished\n";
  return 0;
}

int PliGauss::test(){
  cout << "test mode dayo.\n";
  return 0;
}
int PliGauss::assign_trajtrans(){
  cout << "assign mdcctraj"<<endl;
  cout << "load_gaussian_mixtures()"<<endl;
  Read reader(cfg.fn_interactions);
  int n_frames = reader.loadKKTrajTransHeader(cfg.skip_data);
  int dim = reader.get_dim();
  map<vector<int>, GaussianMixture> gm;
  Read(cfg.fn_gaussians).
    load_gaussian_mixtures(1, cfg.target_columns.size()*dim, gm,
			   cfg.skip_header_gaussian);

  cout << "gm.size():" << gm.size() <<endl;
  cout << "type " << cfg.gmm_type[0] << endl;
  map<vector<int>, GaussianMixture>::iterator itr_gm;
  itr_gm = gm.find(cfg.gmm_type);
  if(itr_gm == gm.end()){
    cerr << "Gaussian mixture was not found: ";
    vector<int>::iterator itr;
    for(itr=cfg.gmm_type.begin();
	itr != cfg.gmm_type.end(); itr++){
      cerr << *itr << " ";
    }
    cerr << endl;
    exit(1);
  }
  vector<Gaussian>::const_iterator itr_gg;
  cout << "Gaussian-ID   average" << endl;
  for(itr_gg = itr_gm->second.get_gauss().begin();
      itr_gg != itr_gm->second.get_gauss().end(); itr_gg++){
    cout << (*itr_gg).get_id();
    for(int i=0; i < dim; i++)
      cout << " " <<  (*itr_gg).get_mu()[i];
    cout << endl;
  }
  
  cout << "loop start" <<endl;
  double **table;
  int n_row = (n_frames - cfg.skip_header) / cfg.skip_data;

  cout << "READ: rows:" << n_row << " cols:" << (int)cfg.target_columns.size()*3 << " fn:" << cfg.fn_interactions<< endl;
  table = new double*[n_row];
  for(int i=0; i < n_row; i++){
    table[i] = new double[(int)cfg.target_columns.size()*3];
  }
  table[0][0] = 1.0;
  //cout << table[0][0] << endl;
  n_row = reader.loadKKTrajTrans(table, cfg.target_columns, n_frames, cfg.skip_data, cfg.skip_header);
  //cout << table[0][0] << endl;

  cout << "nrow:" <<n_row<<endl;

  int n_gauss = (int)itr_gm->second.get_gauss().size();
  int n_col = n_gauss*3;
  cout << "ngauss:" <<n_gauss << " ncol:"<<n_col<<endl;

  double **assign = new double*[n_col];
  for(int i=0; i<n_col; i++){
    assign[i] = new double[n_row];
  }
  stringstream ss;
  //skip header lines
  int i_atom=0;
  vector<int>::iterator itr_a;

  for(itr_a = cfg.target_columns.begin();
      itr_a != cfg.target_columns.end(); itr_a++, i_atom++){

    for(int i_row=0; i_row < n_row; i_row++ ){
      vector<Gaussian>::const_iterator itr_g;
      ublas::vector<double> dat(dim);
      for(int d = 0; d < dim; d++){
	dat[d] = table[i_row][i_atom*dim+d];
      }
      int i_gc=0;
      for(itr_g = itr_gm->second.get_gauss().begin();
	  itr_g != itr_gm->second.get_gauss().end(); itr_g++, i_gc++){
	double dist = (*itr_g).cal_gauss_maharanobis_dist(dat);
	//cout << "dist " << dist << endl;
	double prob_dens = (*itr_g).cal_prob(dat);
	//cout << "dens " << prob_dens << endl;
	double prob_dens_gmm = prob_dens * itr_gm->second.get_pi()[i_gc];
	//cout << "assign : " << i_gc << " - " << i_row <<endl;
	assign[i_gc*3][i_row] = dist;
	assign[i_gc*3+1][i_row] = prob_dens;
	assign[i_gc*3+2][i_row] = prob_dens_gmm;
	//cout << "dist:" << dist << " prob_dens:" << prob_dens << " prob_dens_gmm:" << prob_dens_gmm << endl;
	//cout << "assigned" <<  endl;
      }
    }
  }
  for(int i=0;i<n_row;i++){
    delete[] table[i];
  }
  delete[] table;
  cout << "Output: " << cfg.fn_result<< endl;
  cout<<n_col<<endl;
  Write writer(cfg.fn_result);
  if(cfg.format_output==FMT_ASCII){
    writer.write_assign_data_ascii(assign, n_col, n_row);
  }else if (cfg.format_output==FMT_BIN){
    writer.write_assign_data(assign, n_col, n_row);
  }
  
  //delete variables
  for(int i=0; i<n_col; i++){
    delete[] assign[i];
  }
  delete[] assign;

  return 0;
}

int PliGauss::assign_datatable(){
  cout << "assign_datatable()"<<endl;
  cout << "load_gaussian_mixtures()"<<endl;
  map<vector<int>, GaussianMixture> gm;
  Read(cfg.fn_gaussians).
    load_gaussian_mixtures(1, cfg.target_columns.size()*3, gm,
			   cfg.skip_header_gaussian);

  cout << "gm.size():" << gm.size() <<endl;

  map<vector<int>, GaussianMixture>::iterator itr_gm;
  itr_gm = gm.find(cfg.gmm_type);
  if(itr_gm == gm.end()){
    cerr << "Gaussian mixture was not found: ";
    vector<int>::iterator itr;
    for(itr=cfg.gmm_type.begin();
	itr != cfg.gmm_type.end(); itr++){
      cerr << *itr << " ";
    }
    cerr << endl;
    exit(1);
  }
  int n_gauss = (int)itr_gm->second.get_gauss().size();
  int n_col = n_gauss*3;

  Read reader(cfg.fn_interactions);

  int n_row = 0;
  cout << "Counting the number of rows" << endl;
  vector<double> dat;

  reader.open();
  // Skip header
  for(int i=0; i<cfg.skip_header; i++){
    reader.get_line();
  }
  // Count lines
  int cur_line = -1;
  while(reader.load_datatable_line(dat, cfg.target_columns)==0){
    cur_line++;
    if(cur_line%cfg.skip_data==0) n_row++;
    dat.clear();
  }
  reader.close();

  // Allocate memory
  double **assign = new double*[n_col];
  for(int i=0; i<n_col; i++){
    assign[i] = new double[n_row];
  }
  cout << "The number of rows : " << n_row << endl;
  cout << "Loop start" <<endl;
  reader.open();

  //Write writer(cfg.fn_result);
  //writer.open();

  //skip header lines
  for(int i=0; i<cfg.skip_header; i++){
    reader.get_line();
    //reader.load_datatable_line(dat, cfg.target_columns)
  }

  int i_row = -1;
  cur_line = -1;
  while(reader.load_datatable_line(dat, cfg.target_columns)==0){
    cur_line++;
    if(cur_line%cfg.skip_data!=0){
      dat.clear();
      continue;
    }
    i_row++;
    //vector<double>::iterator itr_tmp;
    //cerr << "aaa ";
    //for(itr_tmp = dat.begin(); itr_tmp != dat.end(); itr_tmp++){
    //cerr << *itr_tmp << " ";
    //}
    //cerr << endl;

    vector<Gaussian>::const_iterator itr_g;
    int i_gc=0;
    for(itr_g = itr_gm->second.get_gauss().begin();
	itr_g != itr_gm->second.get_gauss().end(); itr_g++, i_gc++){
      double dist = (*itr_g).cal_gauss_maharanobis_dist(dat);
      double prob_dens = (*itr_g).cal_prob(dat);
      double prob_dens_gmm = prob_dens * itr_gm->second.get_pi()[i_gc];
      assign[i_gc*3][i_row] = dist;
      assign[i_gc*3+1][i_row] = prob_dens;
      assign[i_gc*3+2][i_row] = prob_dens_gmm;
      //prob_dens_gmm >= cfg.assign_min_prob_dens){
      //ss << dist << "\t" << prob_dens 
      //<< "\t" << prob_dens_gmm << "\t";
      //}
    }
    //ss << endl;
    //writer.write_line(ss.str());
    dat.clear();
  }
  reader.close();
  cout << "Output: " << cfg.fn_result<< endl;
  cout<<n_col<<endl;
  Write writer(cfg.fn_result);
  if(cfg.format_output==FMT_ASCII){
    writer.write_assign_data_ascii(assign, n_col, n_row);
  }else if (cfg.format_output==FMT_BIN){
    writer.write_assign_data(assign, n_col, n_row);
  }
  //delete variables
  for(int i=0; i<n_col; i++){
    delete[] assign[i];
  }
  delete[] assign;
  //writer.close();
  return 0;
}

int PliGauss::assign_tti(){
  cout << "assign_tti()"<<endl;
  
  cout << "load_gauss()"<<endl;
  map<vector<int>, GaussianMixture> gm;
  
  Read(cfg.fn_gaussians).
    load_gaussian_mixtures(6, 3, gm,
			   cfg.skip_header_gaussian);
  
  cout << "gm.size():" << gm.size() <<endl;

  Read reader(cfg.fn_interactions);
  reader.open();
  Write writer(cfg.fn_result);
  writer.open();
  
  TriTriInact tti("",-1,vector<int>(),-1,-1,-1,-1,Coord(),Coord(),Coord());
  int i=-1;
  while(reader.load_tri_tri_interaction_line(tti)==0){
    i++;
    if(cfg.interactions_n_end > 0){
      if(cfg.interactions_n_begin > i) continue;
      else if(cfg.interactions_n_end <= i) break;
    }

    vector<double> x;
    x.push_back(tti.get_ig1().get_x());
    x.push_back(tti.get_ig1().get_y());
    x.push_back(tti.get_ig1().get_z());

    map<vector<int>, GaussianMixture>::iterator itr_gm;
    itr_gm = gm.find(tti.get_type());
    if(itr_gm == gm.end()) continue;
    vector<Gaussian>::const_iterator itr_g;
    int i_gc=0;
    for(itr_g = itr_gm->second.get_gauss().begin();
	itr_g != itr_gm->second.get_gauss().end(); itr_g++, i_gc++){
      double dist = (*itr_g).cal_gauss_maharanobis_dist(x);
      double prob_dens = (*itr_g).cal_prob(x);
      double prob_dens_gmm = prob_dens * itr_gm->second.get_pi()[i_gc];
      if(dist <= cfg.assign_max_distance &&
	 prob_dens_gmm >= cfg.assign_min_prob_dens){
	stringstream ss;
	ss << tti.get_code() <<"\t"
	   << tti.get_id() <<"\t"
	   << (*itr_g).get_id() <<"\t"
	   << dist << "\t" << prob_dens 
	   << "\t" << prob_dens_gmm << endl;
	writer.write_line(ss.str());
      }
    }
  }
  reader.close();
  writer.close();
  return 0;
}
int PliGauss::assign_tai(){
  cout << "assign_tai()"<<endl;
  
  cout << "load_gauss()"<<endl;
  map<vector<int>, GaussianMixture> gm;
  Read(cfg.fn_gaussians).
    load_gaussian_mixtures(4, 3, gm, cfg.skip_header_gaussian);
  cout << "gm.size():" << gm.size() <<endl;

  Read reader(cfg.fn_interactions);
  reader.open();
  Write writer(cfg.fn_result);
  writer.open();
  
  TriAtomInact tai("",-1,vector<int>(),-1,-1,-1,-1,Coord());
  int i=-1;
  while(reader.load_tri_atom_interaction_line(tai, cfg.data_type_cols)==0){
    i++;
    if(cfg.interactions_n_end > 0){
      if(cfg.interactions_n_begin > i) continue;
      else if(cfg.interactions_n_end <= i) break;
    }

    vector<double> x;
    x.push_back(tai.get_ig().get_x());
    x.push_back(tai.get_ig().get_y());
    x.push_back(tai.get_ig().get_z());

    map<vector<int>, GaussianMixture>::iterator itr_gm;
    itr_gm = gm.find(tai.get_type());
    if(itr_gm == gm.end()) continue;
    vector<Gaussian>::const_iterator itr_g;
    int i_gc=0;
    for(itr_g = itr_gm->second.get_gauss().begin();
	itr_g != itr_gm->second.get_gauss().end(); itr_g++, i_gc++){
      double dist = (*itr_g).cal_gauss_maharanobis_dist(x);
      double prob_dens = (*itr_g).cal_prob(x);
      double prob_dens_gmm = prob_dens * itr_gm->second.get_pi()[i_gc];
      if(dist <= cfg.assign_max_distance &&
	 prob_dens_gmm >= cfg.assign_min_prob_dens){
	stringstream ss;
	ss << tai.get_code() << "\t"
	   << tai.get_id() << "\t"
	   << (*itr_g).get_id() <<"\t"
	   << dist << "\t"
	   << prob_dens << "\t" 
	   << prob_dens_gmm << "\t" << endl;
	writer.write_line(ss.str());
      }
    }
  }
  reader.close();
  writer.close();
  return 0;
}

