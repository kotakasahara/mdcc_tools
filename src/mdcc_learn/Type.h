#ifndef __TYPE_H__
#define __TYPE_H__

#include <iostream>
using namespace std;

class TriType{
 public:
  int a1;  int a2;  int a3;
  char b12;  char b23;
  TriType(int inA1,int inA2,int inA3,char inB12,char inB23);
  bool operator==(TriType op);
};
bool operator<(TriType op1,TriType op2);

class ProtriType{
 public:
  int res;  int a1;  int a2;  int a3;
  ProtriType(int inRes,int inA1,int inA2,int inA3);
};
bool operator<(ProtriType op1,ProtriType op2);

#endif 
