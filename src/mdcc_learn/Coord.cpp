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
Coord Coord::operator=(Coord op){
  moveTo(op.getX(), op.getY(), op.getZ());
  return *this;
}
Coord Coord::operator=(Coord *op){
  moveTo(op->getX(), op->getY(), op->getZ());
  return *this;
}
Coord Coord::operator+(Coord op){
  Coord ret;
  ret.moveTo(x+op.getX(), y+op.getY(), z+op.getZ());
  return ret;
}
Coord Coord::operator+(Coord *op){
  Coord ret;
  ret.moveTo(x+op->getX(), y+op->getY(), z+op->getZ());
  return ret;
}
Coord Coord::operator-(Coord op){
  Coord ret;
  ret.moveTo(x-op.getX(), y-op.getY(), z-op.getZ());
  return ret;
}
Coord Coord::operator-(Coord *op){
  Coord ret;
  ret.moveTo(x-op->getX(), y-op->getY(), z-op->getZ());
  return ret;
}
Coord Coord::operator*(double op){
  Coord ret;
  ret.moveTo(x*op, y*op, z*op);
  return ret;
}
Coord Coord::operator+=(Coord op){
  moveTo(x+op.getX(), y+op.getY(), z+op.getZ());
  return *this;
}
Coord Coord::operator+=(Coord *op){
  moveTo(x+op->getX(), y+op->getY(), z+op->getZ());
  return *this;
}
Coord Coord::operator-=(Coord op){
  moveTo(x-op.getX(), y-op.getY(), z-op.getZ());
  return *this;
}
Coord Coord::operator-=(Coord *op){
  moveTo(x-op->getX(), y-op->getY(), z-op->getZ());
  return *this;
}
Coord Coord::operator*=(double op){
  moveTo(x*op, y*op, z*op);
  return *this;
}
bool Coord::operator==(Coord op){
  return (x==op.getX() && y==op.getY() && z==op.getZ());
}
bool Coord::operator!=(Coord op){
  return (x!=op.getX() || y!=op.getY() || z!=op.getZ());
}

double Coord::distTo(Coord inC){
  double dx,dy,dz,r;
  dx = inC.getX()-x;
  dy = inC.getY()-y;
  dz = inC.getZ()-z;
  r=sqrt(dx*dx+dy*dy+dz*dz);
  return r;
}

double Coord::norm(){
  return sqrt(x*x+y*y+z*z);
}
double Coord::innerProduct(Coord inC){
  return x*inC.getX()+y*inC.getY()+z*inC.getZ();
}
Coord Coord::outerProduct(Coord inC){
  return Coord(y*inC.getZ()-z*inC.getY(),
	       z*inC.getX()-x*inC.getZ(),
	       x*inC.getY()-y*inC.getX());
}
double Coord::angle(Coord inC){
  return acos(innerProduct(inC)/(norm()*inC.norm()));
}
double Coord::angle(Coord a, Coord b, Coord c){
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
    cX+=(*itrC).getX();
    cY+=(*itrC).getY();
    cZ+=(*itrC).getZ();
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
    cX+=(*itrC)->getX();
    cY+=(*itrC)->getY();
    cZ+=(*itrC)->getZ();
  }
  cX/=crds.size();
  cY/=crds.size();
  cZ/=crds.size();
  moveTo(cX,cY,cZ);
  return *this;
}
bool operator<(Coord op1,Coord op2){
  if(op1.getX()!=op2.getX())
    return op1.getX()<op2.getX();
  else if(op1.getY()!=op2.getY())
    return op1.getY()<op2.getY();
  return op1.getZ()<op2.getZ();
}

double Coord::calIgTheta(){
  return Coord(0,1,0).angle(*this);
}
double Coord::calIgPhi(){
  Coord crdXZ(x,0,z);
  double phi=Coord(0,0,1).angle(crdXZ);
  if(Coord(0,0,1).outerProduct(crdXZ).getY()>0) phi*=-1;
  return phi;
}
