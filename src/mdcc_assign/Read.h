#ifndef __READ_H__
#define __READ_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <list>
#include <set>
#include <cstdlib>
#include <cstring>
#include "Inact.h"
#include "Gaussian.h"
using namespace std;


class Read{
 private:
  string filename;
  bool op;
  int cur_line;
  bool conv_endian;
  int size_real;
 public:
  ifstream ifs;
  Read(string inFn);
  string getFn(){return filename;};
  bool is_open(){return op;};
  int open();
  int close();
  vector<string> load_config();
  vector<int> loadIntegers();
  vector<string> loadStrings();
  int load_tri_tri_interactions(vector<TriTriInact>& tti,
				int n_begin, int n_end);
  int seek_line(int to_line);
  string get_line();
  int load_tri_tri_interaction_line(TriTriInact& tti);
  int load_tri_atom_interaction_line(TriAtomInact& tai);
  int load_gaussian_mixtures(int n_column_type, int dimension,
			     map<vector<int>, GaussianMixture>& gm,
			     int skip_header_gaussian);
  
  int load_datatable_line(vector<double>& dat,vector<int> column);
  int loadKKTrajTransHeader(int skip);
  int loadKKTrajTrans(double** table, const vector<int> &atomid, int n_frames, int skip, int skip_header);
  bool isConvEndian(){return conv_endian;};
  void setConvEndianTrue(){conv_endian=true;};
  void setConvEndianFalse(){conv_endian=false;};
  template <typename TYPE> int readBinValues(TYPE *recept, int len);
};
inline int reverseEndian(int value){
  int v = value;
  //memcpy(&v, &value, sizeof(v));
  v = (int)((v & 0x00FF00FF) << 8 | (v & 0xFF00FF00) >> 8);
  v = (int)((v & 0x0000FFFF) << 16 | (v & 0xFFFF0000) >> 16);  
  //memcpy(&value, &v, sizeof(value));
  return v;
}

inline float reverseEndian(float value){
  unsigned char v1[4];
  unsigned char v2[4];
  memcpy(&v1, &value, sizeof(v1));

  memcpy(&v2[0],&v1[3], sizeof(v1[3]));
  memcpy(&v2[1],&v1[2], sizeof(v1[2]));
  memcpy(&v2[2],&v1[1], sizeof(v1[1]));
  memcpy(&v2[3],&v1[0], sizeof(v1[0]));

  float ret;
  memcpy(&ret, &v2, sizeof(v2));  
  return ret;
}
inline double reverseEndian(double value){
  unsigned char v[8];
  unsigned char v2[8];
  memcpy(&v, &value, sizeof(v));
  for(int i=0; i < sizeof(v); i++){
    memcpy(&v2[7-i],&v[i], sizeof(v[i]));
  }
  //v = ((v&0x00FF00FF00FF00FF) << 8 | (v & 0xFF00FF00FF00FF00) >> 8);
  //v = ((v & 0x0000FFFF0000FFFF) << 16 | (v & 0xFFFF0000FFFF0000) >> 16);  
  //v = ((v & 0x00000000FFFFFFFF) << 32 | (v & 0xFFFFFFFF00000000) >> 32);  
  double ret;
  memcpy(&ret, v2, sizeof(v2));
  return ret;
}
  

#endif
