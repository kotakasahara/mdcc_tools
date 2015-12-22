#ifndef __FEATURE_DEFINITION_H__
#define __FEATURE_DEFINITION_H__

#include <string>
using namespace std;

enum{
  FTYPE_DUM=0,
  FTYPE_REAL,
  FTYPE_DUM2
};

class FeatureDefinition{
 private:
  int type;
  int column;
  string header;
 public:
  FeatureDefinition(int inType, int inColumn);
  FeatureDefinition(int inType, int inColumn, string inH);
  void setHeader(string inH){header=inH;};
  int getType(){return type;};
  int getColumn(){return column;};
  string getHeader(){return header;};
};

#endif
