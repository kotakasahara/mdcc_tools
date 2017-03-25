#include "Config.h"

Config::Config(){
  mode=M_TEST;
  fn_cfg = "";
  fn_result = "mdccassign_result.txt";
  fn_interactions = "interactions.txt";
  fn_gaussians = "gaussians.txt";
  assign_max_distance = 2.5;
  assign_min_prob_dens = 0.001;
  flg_interaction_type = TTI_TYPE_SYBYL_SYBYL;
  interactions_n_begin = 0;

  skip_data = 1;
  skip_header = 0;
  skip_header_gaussian = 0;

  format_output = FMT_ASCII;
}

Config::~Config(){
}

void Config::set_all(int argn, char* argv[]){
  vector<string> arg;
  int i;
  for(i=1;i<argn;i++)
    arg.push_back(string(argv[i]));
  set_all(arg);
}

void Config::set_all(const vector<string>& arg){
  vector<string>::const_iterator itr;
  string type,val;
  for(itr=arg.begin(); itr!=arg.end(); itr++){
    if(*itr=="-mode"){
      itr++;
      if (*itr == "test")            { mode=M_TEST; }
      else if(*itr == "assign-tti")    { mode=M_ASSIGN_TTI; }
      else if(*itr == "assign-tai")    { mode=M_ASSIGN_TAI; }
      else if(*itr == "assign-table")    { mode=M_ASSIGN_TABLE; }
      else if(*itr == "assign-mdcctraj")    { mode=M_ASSIGN_TRAJTRANS; }
      else{
	cerr<<"invalid mode [" << *itr << "]\n"; exit(1);
      }
    }
    else if(*itr=="-fn-cfg"){ fn_cfg=*++itr; }
    else if(*itr=="-fn-interactions" || *itr=="-fn-data-table"){ fn_interactions=*++itr; }
    else if(*itr=="-fn-gaussians"){ fn_gaussians=*++itr; }
    else if(*itr=="-fn-result"){ fn_result=*++itr; }
    else if(*itr=="-format-output"){ 
      itr++;
      if (*itr == "binary") {format_output=FMT_BIN;}
      else if(*itr == "ascii") { format_output=FMT_ASCII;}
      else{
	cerr << "invalid format-output: " << *itr << endl;
	exit(1);
      }
    }
    else if(*itr=="-assign-max-dist"){ assign_max_distance=atof((*++itr).c_str()); }
    else if(*itr=="-assign-min-prob-dens"){ assign_min_prob_dens=atof((*++itr).c_str()); }
    else if(*itr=="-interactions-range"){ 
      interactions_n_begin=atoi((*++itr).c_str());
      interactions_n_end=atoi((*++itr).c_str());
    }else if(*itr=="-target-column"){ target_columns.push_back(atoi((*++itr).c_str())); }
    else if(*itr=="-gmm-type"){ gmm_type.push_back(atoi((*++itr).c_str()));}
    else if(*itr=="-skip-data"){ skip_data=atoi((*++itr).c_str()); }
    else if(*itr=="-skip-header"){ skip_header=atoi((*++itr).c_str()); }
    else if(*itr=="-skip-header-gaussian"){ skip_header_gaussian=atoi((*++itr).c_str()); }
    else if(*itr=="-data-type-col"){ data_type_cols.push_back(atoi((*++itr).c_str()));}
    else{
      cerr << "unknown keyword <" << *itr << ">" << endl;
    }
  }
}

