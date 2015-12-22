#ifndef __TYPE_H__
#define __TYPE_H__

#include <iostream>
using namespace std;
#include "define.h"

class TriType{
 public:
  int a1;  int a2;  int a3;
  char b12;  char b23;
  TriType();
  TriType(int inA1,int inA2,int inA3,char inB12,char inB23);
  TriType(int inA1,int inA2,int inA3);
  string getTypeSt();
  bool operator==(TriType op);
};
bool operator<(TriType op1,TriType op2);

class ProtriType{
 public:
  int flgType; // 0:atom names,  1:sybyl types
  int res;  int a1;  int a2;  int a3;
  ProtriType();
  ProtriType(int inRes,int inA1,int inA2,int inA3,int flgPtriType);
  ProtriType(int inRes,int inA1,int inA2,int inA3);
  string getTypeSt();
  string getAtomNameSt();
  string getSybylTypeSt();
  bool operator==(ProtriType op);
};
bool operator<(ProtriType op1,ProtriType op2);

#endif 
