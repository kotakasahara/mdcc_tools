#include "Coord.h"
using namespace std;

Coord::Coord(){
}
Coord::Coord(double inX, double inY, double inZ){
  moveTo(inX, inY, inZ);
}
Coord::Coord(vector<Coord> crds){
  moveToCenter(crds);
}
Coord::Coord(vector<Coord*> crds){
  moveToCenter(crds);
}
Coord Coord::operator=(const Coord& op) {
  moveTo(op.get_x(), op.get_y(), op.get_z());
  return *this;
}
//Coord Coord::operator=(const Coord& *op) const{
//  moveTo(op->get_x(), op->get_y(), op->get_z());
//  return *this;
//}
Coord Coord::operator+(const Coord& op) const {
  Coord ret;
  ret.moveTo(x+op.get_x(), y+op.get_y(), z+op.get_z());
  return ret;
}
//Coord Coord::operator+(const Coord& *op) const{
//  Coord ret;
//  ret.moveTo(x+op->get_x(), y+op->get_y(), z+op->get_z());
//  return ret;
//}
Coord Coord::operator-(const Coord& op) const {
  Coord ret;
  ret.moveTo(x-op.get_x(), y-op.get_y(), z-op.get_z());
  return ret;
}
//Coord Coord::operator-(const Coord& *op) const{
//  Coord ret;
//  ret.moveTo(x-op->get_x(), y-op->get_y(), z-op->get_z());
//  return ret;
//}
Coord Coord::operator*(const double& op) const {
  Coord ret;
  ret.moveTo(x*op, y*op, z*op);
  return ret;
}
Coord Coord::operator+=(const Coord& op){
  moveTo(x+op.get_x(), y+op.get_y(), z+op.get_z());
  return *this;
}
//Coord Coord::operator+=(const Coord& *op) const {
//  moveTo(x+op->get_x(), y+op->get_y(), z+op->get_z());
//  return *this;
//}
Coord Coord::operator-=(const Coord& op){
  moveTo(x-op.get_x(), y-op.get_y(), z-op.get_z());
  return *this;
}
//Coord Coord::operator-=(const Coord& *op) const{
//  moveTo(x-op->get_x(), y-op->get_y(), z-op->get_z());
//  return *this;
//}
Coord Coord::operator*=(double op){
  moveTo(x*op, y*op, z*op);
  return *this;
}
bool Coord::operator==(const Coord& op) const{
  return (x==op.get_x() && y==op.get_y() && z==op.get_z());
}
bool Coord::operator!=(const Coord& op) const{
  return (x!=op.get_x() || y!=op.get_y() || z!=op.get_z());
}

double Coord::distTo(const Coord& inC) const {
  double dx,dy,dz,r;
  dx = inC.get_x()-x;
  dy = inC.get_y()-y;
  dz = inC.get_z()-z;
  r=sqrt(dx*dx+dy*dy+dz*dz);
  return r;
}

double Coord::norm() const{
  return sqrt(x*x+y*y+z*z);
}
double Coord::innerProduct(const Coord& inC) const {
  return x*inC.get_x()+y*inC.get_y()+z*inC.get_z();
}
Coord Coord::outerProduct(const Coord& inC) const {
  return Coord(y*inC.get_z()-z*inC.get_y(),
	       z*inC.get_x()-x*inC.get_z(),
	       x*inC.get_y()-y*inC.get_x());
}
double Coord::angle(const Coord& inC) const {
  return acos(innerProduct(inC)/(norm()*inC.norm()));
}
double Coord::angle(const Coord& a, const Coord& b, const Coord& c) const{
  Coord AB = b-a;
  Coord BC = c-b;
  return acos(AB.innerProduct(BC)/(AB.norm()*BC.norm()));
}
double Coord::torsion(Coord a, Coord b, Coord c, Coord d){
  Coord AB = b-a;
  Coord BC = c-b;
  Coord CD = d-c;
  Coord opABC = AB.outerProduct(BC);
  Coord opBCD = BC.outerProduct(CD);
  return opABC.angle(opBCD);
}
void Coord::rotateX(double theta){
  double tmpY=y;
  y=y*cos(theta)-z*sin(theta);
  z=tmpY*sin(theta)+z*cos(theta);
}
void Coord::rotateY(double theta){
  double tmpX=x;
  x=x*cos(theta)+z*sin(theta);
  z=-tmpX*sin(theta)+z*cos(theta);
}
void Coord::rotateZ(double theta){
  double tmpX=x;
  x=x*cos(theta)-y*sin(theta);
  y=tmpX*sin(theta)+y*cos(theta);
}
Coord Coord::moveToCenter(vector<Coord> crds){
  double cX=0;
  double cY=0;
  double cZ=0;
  vector<Coord>::iterator itrC;
  for(itrC=crds.begin();itrC!=crds.end();itrC++){
    cX+=(*itrC).get_x();
    cY+=(*itrC).get_y();
    cZ+=(*itrC).get_z();
  }
  cX/=crds.size();
  cY/=crds.size();
  cZ/=crds.size();
  moveTo(cX,cY,cZ);
  return *this;
}
Coord Coord::moveToCenter(vector<Coord*> crds){
  double cX=0;
  double cY=0;
  double cZ=0;
  vector<Coord*>::iterator itrC;
  for(itrC=crds.begin();itrC!=crds.end();itrC++){
    cX+=(*itrC)->get_x();
    cY+=(*itrC)->get_y();
    cZ+=(*itrC)->get_z();
  }
  cX/=crds.size();
  cY/=crds.size();
  cZ/=crds.size();
  moveTo(cX,cY,cZ);
  return *this;
}
bool operator<(Coord op1,Coord op2){
  if(op1.get_x()!=op2.get_x())
    return op1.get_x()<op2.get_x();
  else if(op1.get_y()!=op2.get_y())
    return op1.get_y()<op2.get_y();
  return op1.get_z()<op2.get_z();
}

double Coord::calIgTheta(){
  return Coord(0,0,1).angle(*this);
}
double Coord::calIgPhi(){
  Coord crdXY(x,y,0);
  double phi=Coord(1,0,0).angle(crdXY);
  if(Coord(1,0,0).outerProduct(crdXY).get_z()>0) phi*=-1;
  return phi;
}
