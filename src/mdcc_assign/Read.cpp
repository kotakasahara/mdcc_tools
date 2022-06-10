#include "Read.h"

Read::Read(string inFn){
  op=false;
  filename = inFn;
}

int Read::open(){
  ifs.open(filename.c_str());
  if(!ifs){
    cerr <<"Cannot open "<< filename << "." <<endl;
    return 1;
  }
  op=true;
  return 0;
}
int Read::close(){
  ifs.close();
  op=false;
  return 0;
}

vector<string> Read::load_config(){
  vector<string> vconf;
  open();
  string buf;
  while(ifs && getline(ifs, buf)){
    int pos1 = buf.find_first_of("#;");
    if(pos1 != string::npos){
      buf = buf.substr(0, pos1);
    }
    if(buf.size()==0) continue;
    stringstream ss(buf);
    while(ss>>buf){
      vconf.push_back(buf);
    }
  }

  close();
  return vconf;

}

vector<int> Read::loadIntegers(){
  vector<int> intvec;
  open();
  string buf;
  while(getline(ifs,buf)){
    intvec.push_back(atoi(buf.c_str()));
  }
  close();
  return intvec;
}
vector<string> Read::loadStrings(){
  vector<string> strvec;
  open();
  string buf;
  while(getline(ifs,buf)){
    stringstream ss(buf);
    string str;
    ss>>str;
    strvec.push_back(str);
  }
  close();
  return strvec;
}
int Read::seek_line(int to_line){
  string buf;
  while(cur_line < to_line){
    getline(ifs,buf);
    cur_line++;
  }
  return 0;
}
string Read::get_line(){
  string buf;
  getline(ifs,buf);
  return buf;
}


int Read::load_datatable_line(vector<double>& dat,
			      vector<int> column){
  string buf;
  string tmp;

  if(!getline(ifs,buf)) return 1;
  stringstream ss(buf);
  //cout << "TEST * " << ss.str() << endl;
  double tmp_dat;
  int i = 0;
  while(ss>>tmp_dat){
    vector<int>::iterator itr_find;
    for(itr_find=column.begin(); itr_find!=column.end(); itr_find++){
      if(*itr_find == i){
	dat.push_back(tmp_dat);
      }
    }
    //itr_find = find(column.begin(), column.end(), i);
    
    //if(itr_find != column.end()){
    //  dat.push_back(tmp_dat);
    //}
    i++;
  }
  return 0;
}

int Read::load_tri_tri_interaction_line(TriTriInact& tti){
  string buf;

  if(!getline(ifs,buf)){
    return 1;
  }
  
  string code;
  int id;
  int ptriid;
  int resid;
  int restype;
  int ptri1; int ptri2; int ptri3;
  int ltriid;
  int ltri1; int ltri2; int ltri3;
  
  double ig1x,ig1y,ig1z;
  double ig2x,ig2y,ig2z;
  double ig3x,ig3y,ig3z;
  //    double igCx,igCy,igCz;
  //    double igNx,igNy,igNz;
  stringstream ss(buf);
  ss>>code>>id>>ptriid>>resid>>restype
    >>ptri1>>ptri2>>ptri3
    >>ltriid
    >>ltri1>>ltri2>>ltri3
    >>ig1x>>ig1y>>ig1z
    >>ig2x>>ig2y>>ig2z
    >>ig3x>>ig3y>>ig3z;
  //      >>igCx>>igCy>>igCz
  //      >>igNx>>igNy>>igNz;
  vector<int> type;
  type.push_back(ptri1);
  type.push_back(ptri2);
  type.push_back(ptri3);
  type.push_back(ltri1);
  type.push_back(ltri2);
  type.push_back(ltri3);
  TriTriInact tmptti(code, id, type, ptriid, resid, restype,
		     ltriid,
		     Coord(ig1x,ig1y,ig1z),
		     Coord(ig2x,ig2y,ig2z),
		     Coord(ig3x,ig3y,ig3z));
  tti = tmptti;
  return 0;
}
int Read::load_tri_atom_interaction_line(TriAtomInact& tai,
					 vector<int> type_cols){
  string buf;
  if(type_cols.empty()){
    type_cols.push_back(5);
    type_cols.push_back(6);
    type_cols.push_back(7);
    type_cols.push_back(9);
  }

  if(!getline(ifs,buf)){
    return 1;
  }
  
  string code;
  int id;
  int ptriid;
  int resid;
  int restype;
  int ptri1; int ptri2; int ptri3;
  int latomidv;
  int latomtype;
  int ltri1; int ltri2; int ltri3;
  
  double igx,igy,igz;

  stringstream ss(buf);
  string token;
  vector<string> terms;
  while ( ss >> token ){
    terms.push_back(token);
  }
  code = terms[0];
  id = atoi(terms[1].c_str());
  ptriid = atoi(terms[2].c_str());
  resid = atoi(terms[3].c_str());
  restype = atoi(terms[4].c_str());
  latomidv = atoi(terms[8].c_str());
  igx = atof(terms[10].c_str());
  igy = atof(terms[11].c_str());
  igz = atof(terms[12].c_str());
  ptri1 = atoi(terms[type_cols[0]].c_str());
  ptri2 = atoi(terms[type_cols[1]].c_str());
  ptri3 = atoi(terms[type_cols[2]].c_str());
  latomtype  = atoi(terms[type_cols[3]].c_str());
  //cout << "dbg0325 " << ptri1 << "-" << ptri2 << "-" << ptri3 <<"-"<<latomtype <<endl;
  //ss >> code >> id >> ptriid >> resid >> restype
  //>> ptri1 >> ptri2 >> ptri3
  //>>latomidv
  //>>latomtype
  //>>igx>>igy>>igz;

  // >>igCx>>igCy>>igCz
  // >>igNx>>igNy>>igNz;
  vector<int> type;
  type.push_back(ptri1);
  type.push_back(ptri2);
  type.push_back(ptri3);
  type.push_back(latomtype);

  TriAtomInact tmptai(code, id, type, ptriid, resid, restype,
		     latomidv,
		     Coord(igx,igy,igz));
  tai = tmptai;
  return 0;
}
int Read::load_res_atom_interaction_line(ResAtomInact& rai,
					 vector<int> type_cols){
  string buf;
  if(type_cols.empty()){
    type_cols.push_back(3);
    type_cols.push_back(6);
    type_cols.push_back(10);
  }

  if(!getline(ifs,buf)){
    return 1;
  }
  
  string code;
  int id;
  int ptriid;
  int resid;
  int restype;
  int latomidv;
  int bbsc;
  int latomtype;
  
  double igx,igy,igz;

  stringstream ss(buf);
  string token;
  vector<string> terms;
  while ( ss >> token ){
    terms.push_back(token);
  }
  code = terms[0];
  id = atoi(terms[1].c_str());
  resid = atoi(terms[2].c_str());
  restype = atoi(terms[3].c_str());
  latomidv = atoi(terms[4].c_str());
  igx = atof(terms[7].c_str());
  igy = atof(terms[8].c_str());
  igz = atof(terms[9].c_str());
  restype = get_res_id(terms[type_cols[0]]);
  bbsc = atoi(terms[type_cols[1]].c_str());
  latomtype = atoi(terms[type_cols[2]].c_str());
  vector<int> type;
  type.push_back(restype);
  type.push_back(bbsc);
  type.push_back(latomtype);

  ResAtomInact tmprai(code, id, type, resid, restype,
		      latomidv,
		      Coord(igx,igy,igz));
  rai = tmprai;
  return 0;
}
int Read::load_tri_tri_interactions(vector<TriTriInact>& tti,
				    int n_begin,
				    int n_end){
  string buf;
  open();
  int n_line = -1;
  while(getline(ifs,buf)){
    n_line++;
    if(n_end > 0){
      if(n_line < n_begin)     continue;
      else if(n_line >= n_end) break;
    }  

    string code;
    int id;
    int ptriid;
    int resid;
    int restype;
    int ptri1; int ptri2; int ptri3;
    int ltriid;
    int ltri1; int ltri2; int ltri3;

    double ig1x,ig1y,ig1z;
    double ig2x,ig2y,ig2z;
    double ig3x,ig3y,ig3z;
    //    double igCx,igCy,igCz;
    //    double igNx,igNy,igNz;
    stringstream ss(buf);
    ss>>code>>id>>ptriid>>resid>>restype
      >>ptri1>>ptri2>>ptri3
      >>ltriid
      >>ltri1>>ltri2>>ltri3
      >>ig1x>>ig1y>>ig1z
      >>ig2x>>ig2y>>ig2z
      >>ig3x>>ig3y>>ig3z;
      //      >>igCx>>igCy>>igCz
      //      >>igNx>>igNy>>igNz;
    vector<int> type;
    type.push_back(ptri1);
    type.push_back(ptri2);
    type.push_back(ptri3);
    type.push_back(ltri1);
    type.push_back(ltri2);
    type.push_back(ltri3);
    TriTriInact tmptti(code, id, type, ptriid, resid, restype,
		       ltriid,
		       Coord(ig1x,ig1y,ig1z),
		       Coord(ig2x,ig2y,ig2z),
		       Coord(ig3x,ig3y,ig3z));

    tti.push_back(tmptti);
  }
  close();
  return 0;
}

int Read::load_gaussian_mixtures(int n_column_type, int dimension,
				 map<vector<int>, GaussianMixture>& gm,
				 int skip_header_gaussian){
  open();
  vector<int> prev_type;
  GaussianMixture new_gm;
  string buf;

  // skip the header line
  for(int i=0; i<skip_header_gaussian; i++){
    if(!getline(ifs,buf)) return 1;
  } 
  while(getline(ifs,buf)){
    vector<int> type;
    stringstream ss(buf);
    int id;
    string tmp;
    ss >> id;
    for(int i = 0; i<n_column_type; i++){
      ss >> tmp;
      type.push_back(atoi(tmp.c_str()));
    }
    double pi;
    ss >> pi;
    ublas::vector<double> mu(dimension);
    for(int i=0; i<dimension; i++){
      double tmp_mu;
      ss >> tmp_mu;
      mu[i] = tmp_mu;
    }
    ublas::matrix<double> sigma(dimension,dimension);
    for(int i=0; i<dimension; i++){
      for(int j=0; j<dimension; j++){
	double tmp_sigma;
	ss >> tmp_sigma;
	sigma(i,j) = tmp_sigma;
      }
    }
    Gaussian new_ge(id, type, mu, sigma);
    if(!prev_type.empty() && type!=prev_type){
      //      cout<<"read "<<prev_type<<endl;
      gm.insert(make_pair(prev_type,new_gm));
      new_gm = GaussianMixture();
    }
    new_gm.push_mixture_element(pi,new_ge);
    prev_type = type;
  }
  gm.insert(make_pair(prev_type,new_gm));
  close();
  return 0;
}
int Read::loadKKTrajTransHeader(int skip){
  open();
  int magic;
  int n_atoms;
  int n_frames;

  ifs.read((char*)&magic, sizeof(int));
  if(magic != 1993){
    setConvEndianTrue();
    //cerr << magic << endl;
    magic = reverseEndian(magic);
    if(magic != 1993){
      cerr << magic << endl;
      cerr << "ERROR: the first 4 bytes do not indicate 1993" << endl;
      exit(0);
    }
  }

  readBinValues(&size_real, 1);
  //cout <<  size_real << endl;
  readBinValues(&n_atoms, 1);
  //cout << n_atoms << endl;
  readBinValues(&n_frames, 1);
  //cout << n_frames << endl;
  readBinValues(&dim, 1);
  int n_row = n_frames/skip;
  return n_frames;
  
}
int Read::loadKKTrajTrans(double** table, const vector<int> &atomid, int n_frames, int skip, int skip_header){  
  int i_atom = 0;
  int n_row = (n_frames - skip_header)/skip;
  int size_header = 20;
  vector<int>::const_iterator itr_a;
  for(itr_a = atomid.begin(); itr_a != atomid.end(); itr_a++){
    ifs.seekg(size_header);
    for(int i=0; i<(*itr_a)*dim; i++){
      ifs.seekg(size_real * n_frames, ios::cur);
    }
    if(ifs.tellg() < size_header){
      cout << "ERROR?? ifs.tellg()<0"<<endl;
      cerr << "ERROR?? ifs.tellg()<0"<<endl;
      exit(1);
    }
    // ifs.seekg(16 + size_real * n_frames * (*itr_a)*3);
    //cout << "size_real: "<< size_real << " n_frams:" << n_frames <<" " << (*itr_a)*3<<endl;
    cout << "atomid:" << (*itr_a) << " " << ifs.tellg() << " " << (unsigned long long)(size_header + size_real * n_frames * (*itr_a)*3) << endl;
    for(int i_dim=0; i_dim < dim; i_dim++){
      //cout << "dim " << i_dim << endl;
      for(int i_frm=0; i_frm < n_frames; i_frm++){
	int i_row = (i_frm - skip_header)/skip;
	double buf;
	if(size_real==4){
	  float buff;
	  readBinValues(&buff, 1);  
	  //cout << i_frm << " * " << i_dim << " : " << buff << endl;;
	  buf = (double)buff;
	}else if(size_real==8){
	  readBinValues(&buf, 1);  
	}
	
	if(i_frm%skip==0 && i_frm >= skip_header){
	  table[i_row][i_atom*3+i_dim] = buf;
	  //cout << ifs.tellg() << " " << buf*10.0 << endl;
	}
      }
    }
    i_atom++;
  }
  close();
  return n_row;
}

template <typename TYPE> int Read::readBinValues(TYPE *recept, int len){
  ifs.read((char*)recept, sizeof(TYPE)*len);
  if(isConvEndian()){
    for(int i=0;i<len;i++){
      recept[i] = reverseEndian(recept[i]);
    }
  }
  return 0;
}
