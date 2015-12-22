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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "FeatureDefinition.h"
#include "Gaussian.h"
using namespace std;
namespace ublas=boost::numeric::ublas;

class Read{
 private:
  string filename;
  bool op;
  bool conv_endian;
 public:
  ifstream ifs;
  Read(string inFn);
  string getFn(){return filename;};
  bool is_open(){return op;};
  int open();
  int close();
  vector<string> loadConfig();
  vector<int> loadIntegers();
  vector<string> loadStrings();
  ublas::vector<ublas::vector<double> > loadKKTrajTrans(vector<FeatureDefinition> fDef, int skip=1, int skip_header=0, int end_frame=-1);
  ublas::vector<ublas::vector<double> > loadDataTable(vector<FeatureDefinition> fDef, int skip=1, int skip_header=0, int end_frame=-1);
  int loadGaussianMixtures(int n_column_type, int dimension,
			   map<string,GaussianMixture> *gm);
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
  int v;
  memcpy(&v, &value, sizeof(v));
  v = ((v & 0x00FF00FF) << 8 | (v & 0xFF00FF00) >> 8);
  v = ((v & 0x0000FFFF) << 16 | (v & 0xFFFF0000) >> 16);  
  float ret;
  memcpy(&ret, &v, sizeof(ret));
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
