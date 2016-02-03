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
  setConvEndianFalse();
  return 0;
}
int Read::close(){
  ifs.close();
  op=false;
  return 0;
}

vector<string> Read::loadConfig(){
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
  int tmp;
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


ublas::vector<ublas::vector<double> > Read::loadKKTrajTrans(vector<FeatureDefinition> fDef, int skip, int skip_header, int end_frame){
  open();
  int magic;
  int size_real;
  int n_atoms;
  int n_frames;
  int dim;

  ifs.read((char*)&magic, sizeof(int));
  if(magic != 1993){
    setConvEndianTrue();
    ///cerr << magic << endl;
    magic = reverseEndian(magic);
    if(magic != 1993){
      //cerr << magic << endl;
      cerr << "ERROR: the first 4 bytes do not indicate 1993" << endl;
      exit(0);
    }
  }

  readBinValues(&size_real, 1);
  cout <<  "Size of real value: " << size_real << endl;
  readBinValues(&n_atoms, 1);
  cout << "The number of elements : " << n_atoms << endl;
  readBinValues(&n_frames, 1);
  cout << "The number of samples : " << n_frames << endl;
  readBinValues(&dim, 1);
  cout << "Dimension : " << dim << endl;
  int size_header = 20;

  if(end_frame == -1) end_frame = n_frames;
  
  int n_row = (end_frame - skip_header) / skip;
  ublas::vector<ublas::vector<double> > dataTable(n_row);

  for(int i_row=0;i_row<n_row;i_row++){    
    dataTable[i_row]=ublas::vector<double>((int)fDef.size() * dim);
  }
  
  vector<FeatureDefinition>::iterator itrFDef;    
  int i_atom = 0;
  
  for(itrFDef=fDef.begin(); itrFDef!=fDef.end(); itrFDef++){
    //cout << "col_id: " << (*itrFDef).getColumn() << " size_real:" << size_real << " n_frames:" << n_frames <<endl;

    ifs.seekg(size_header); // size of the header
    for(int i=0; i<(*itrFDef).getColumn() * dim; i++){
      ifs.seekg(size_real * n_frames, ios::cur);
    }

    //cout << (*itrFDef).getColumn() << endl;
    for(int i_dim=0; i_dim < dim; i_dim++){
      //cout << "dim " << i_dim << endl;
      for(int i_frm=0; i_frm<n_frames; i_frm++){
	double buf;
	if(size_real==4){
	  float buff;
	  readBinValues(&buff, 1);  
	  buf = (double)buff;
	}else if(size_real==8){
	  readBinValues(&buf, 1);  
	}
	if(i_frm%skip==0 && i_frm >= skip_header && i_frm <= end_frame){
	  int i_row = (i_frm - skip_header)/skip;
	  
	  dataTable[i_row][i_atom*dim+i_dim] = buf;
	  //cout << "d " << i_frm << " " << i_row << " " << buf << endl;
	  //cout << ifs.tellg() << " " << buf << endl;
	}
      }
    }
    i_atom++;
  }  
  close();
  return dataTable;
}

ublas::vector<ublas::vector<double> > Read::loadDataTable(vector<FeatureDefinition> fDef, int skip, int skip_header, int end_frame){
  vector<vector<double> > table;
  //  cout << "dbg-1"<<endl;
  open();
  //  cout << "dbg-0.5"<<endl;
  string buf;
  int nCol=-1;
  int nLine=0;
  //cout << "dbg0 "<< skip_header <<endl;
  
  for(int i=0;i<skip_header;i++){
    getline(ifs,buf);
    //cout<<"skip!" <<endl;
  }
  while(getline(ifs,buf) && (end_frame==-1 || nLine <= end_frame)){
    nLine++;
    if (nLine%skip != 0) continue;
    stringstream ss(buf);
    if(buf[0] == '#') continue;
    vector<string> tmpSt;
    string str;
    while(ss>>str){
      tmpSt.push_back(str);
    }
    vector<double> row;
    vector<FeatureDefinition>::iterator itrFDef;
    for(itrFDef=fDef.begin(); itrFDef!=fDef.end(); itrFDef++){
      row.push_back(atof(tmpSt[(*itrFDef).getColumn()].c_str()));
    }
    if(nCol==-1) nCol=row.size();
    else if(nCol>(int)row.size()){
      cerr<<"Invalid value, line:"<<nLine<<endl;
      return ublas::vector<ublas::vector<double> >();
    }
    table.push_back(row);
    //    cout << "nline "<<nLine<<endl;
  }
  cout << "nrow:"<<table.size()<<" ncol:"<<nCol<<endl;
  ublas::vector<ublas::vector<double> > dataTable(table.size());
  for(int i=0;i<(int)dataTable.size();i++){
    dataTable[i]=ublas::vector<double>(nCol);
    for(int j=0;j<nCol;j++){
      dataTable[i][j]=table[i][j];
      //cout << dataTable[i][j] << endl;
    }
  }
  close();
  return dataTable;
}

int Read::loadGaussianMixtures(int n_column_type, int dimension,
			       map<string,GaussianMixture> *gm){
  open();
  string prev_type="";
  GaussianMixture new_gm;
  string buf;
  while(getline(ifs,buf)){
    stringstream type;
    stringstream ss(buf);
    int id;
    string tmp;
    ss >> id;
    ss >> tmp;
    type << tmp;
    for(int i = 1; i<n_column_type; i++){
      ss >> tmp;
      type << "\t" << tmp;
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
    Gaussian new_ge(id, type.str(), mu, sigma);
    if(prev_type!="" && type.str()!=prev_type){
      //      cout<<"read "<<prev_type<<endl;
      gm->insert(make_pair(prev_type,new_gm));
      new_gm = GaussianMixture();
    }
    new_gm.pushMixtureElement(pi,new_ge);
    prev_type = type.str();
  }
  gm->insert(make_pair(prev_type,new_gm));
  close();
  return 0;
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
