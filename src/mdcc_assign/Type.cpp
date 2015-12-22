#include "Type.h"
TriType::TriType(){
}
TriType::TriType(int inA1,int inA2,int inA3,char inB12,char inB23){
  a1=inA1; a2=inA2; a3=inA3; b12=inB12; b23=inB23;
}
TriType::TriType(int inA1,int inA2,int inA3){
  a1=inA1; a2=inA2; a3=inA3;
}
string TriType::getTypeSt(){
  return sybNew[a1]+"-"+sybNew[a2]+"-"+sybNew[a3];
}
bool TriType::operator==(TriType op){
  return a1==op.a1 &&
    a2==op.a2 && a3==op.a3;
}
bool operator<(TriType op1,TriType op2){
  if(op1.a1!=op2.a1) return op1.a1<op2.a1;
  else if(op1.a2!=op2.a2) return op1.a2<op2.a2;
  else if(op1.a3!=op2.a3) return op1.a3<op2.a3;
  else if(op1.b12!=op2.b12) return op1.b12<op2.b12;
  else return op1.b23<op2.b23;
}
ProtriType::ProtriType(){
}
ProtriType::ProtriType(int inRes,int inA1,int inA2,int inA3, int flgPtriType){
  res=inRes; a1=inA1; a2=inA2; a3=inA3;
  flgType = flgPtriType;
}
ProtriType::ProtriType(int inRes,int inA1,int inA2,int inA3){
  res=inRes; a1=inA1; a2=inA2; a3=inA3;
}
string ProtriType::getTypeSt(){
  if(flgType==1) return getSybylTypeSt();
  return getAtomNameSt();
}
string ProtriType::getAtomNameSt(){
  return aminoacid[res]+
    "-"+aaatomname_serial[a1]+
    "-"+aaatomname_serial[a2]+
    "-"+aaatomname_serial[a3];
}
string ProtriType::getSybylTypeSt(){
  return sybNew[a1]+
    "-"+sybNew[a2]+
    "-"+sybNew[a3];
}
bool operator<(ProtriType op1,ProtriType op2){
  if(op1.res!=op2.res) return op1.res<op2.res;
  else if(op1.a1!=op2.a1) return op1.a1<op2.a1;
  else if(op1.a2!=op2.a2) return op1.a2<op2.a2;
  else return op1.a3<op2.a3;
}
bool ProtriType::operator==(ProtriType op){
  return res==op.res && a1==op.a1 &&
    a2==op.a2 && a3==op.a3;
}
