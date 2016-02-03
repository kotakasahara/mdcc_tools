#include "Write.h"

Write::Write(string inFn){
  op=false;
  filename =inFn;
}
int Write::open(){
  ofs.open(filename.c_str());
  if(!ofs){
    cerr <<"Cannot open "<<filename<<"."<<endl;
    return 1;
  }
  op=true;
  return 0;
}
int Write::openApp(){
  ofs.open(filename.c_str(),ios::app);
  if(!ofs){
    cerr <<"Cannot open "<<filename<<"."<<endl;
    return 1;
  }
  op=true;
  return 0;
}
int Write::close(){
  ofs.close();
  op=false;
  return 0;
}

int Write::write_assigned_tti(){
  open();
  close();
  return 0;
}

int Write::write_line(const string& st){
  ofs<<st;
  return 0;
}


int Write::write_assign_data(double** assign, int n_col, int n_row){
  open();
  int magic=1993;
  ofs.write((const char*)&magic, sizeof magic);
  ofs.write((const char*)&n_col, sizeof n_col);
  ofs.write((const char*)&n_row, sizeof n_row);
  for(int i_col=0; i_col < n_col; i_col++){
    for(int i_row=0; i_row < n_row; i_row++){
      ofs.write((const char*)&assign[i_col][i_row], sizeof(double));
    }
  }
  close();
  return 0;
}
int Write::write_assign_data_ascii(double** assign, int n_col, int n_row){
  open();

  for(int i_row=0; i_row < n_row; i_row++){
    stringstream ss;
    for(int i_col=0; i_col < n_col; i_col++){
      ss << assign[i_col][i_row];
      if(i_col != n_col-1) ss << "\t";
    }
    ss << endl;
    ofs << ss.str();
  }
  close();
  return 0;
}
