#ifndef __COORD_H__
#define __COORD_H__

#include <cmath>
#include <vector>
#include "define.h"
using namespace std;

class Coord {
 private:
  double x,y,z;
 public:
  Coord();
  Coord(double inX, double inY, double inZ);
  Coord(vector<Coord> crds);
  Coord(vector<Coord*> crds);
  Coord operator=(Coord op);
  Coord operator=(Coord *op);
  Coord operator+(Coord op);
  Coord operator+(Coord *op);
  Coord operator-(Coord op);
  Coord operator-(Coord *op);
  Coord operator*(double op);
  Coord operator+=(Coord op);
  Coord operator+=(Coord *op);
  Coord operator-=(Coord op);
  Coord operator-=(Coord *op);
  Coord operator*=(double op);
  bool operator==(Coord op);
  bool operator!=(Coord op);
  double getX() {return x;};
  double getY() {return y;};
  double getZ() {return z;};
  void setX(double inX){x=inX;};
  void setY(double inY){y=inY;};
  void setZ(double inZ){z=inZ;};
  void moveTo(double inX, double inY, double inZ){x=inX; y=inY; z=inZ;};
  void moveTo(Coord inC){x=inC.getX(); y=inC.getY(); z=inC.getZ();};
  double distTo(Coord inC);
  double norm();
  double innerProduct(Coord inC);
  Coord outerProduct(Coord inC);
  double angle(Coord a, Coord b, Coord c);
  double angle(Coord inC);
  double torsion(Coord a, Coord b, Coord c, Coord d);
  void rotateX(double theta);
  void rotateY(double theta);
  void rotateZ(double theta);
  Coord moveToCenter(vector<Coord> crds);
  Coord moveToCenter(vector<Coord*> crds);
  double calIgTheta();
  double calIgPhi();
};
bool operator<(Coord op1,Coord op2);
#endif
