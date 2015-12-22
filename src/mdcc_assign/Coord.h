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
  Coord operator=(const Coord& op);
  //  Coord operator=(const Coord& *op) const;
  Coord operator+(const Coord& op) const;
  //  Coord operator+(const Coord& *op) const;
  Coord operator-(const Coord& op) const;
  //  Coord operator-(const Coord& *op) const;
  Coord operator*(const double& op) const;
  Coord operator+=(const Coord& op);
  //  Coord operator+=(const Coord& *op) const;
  Coord operator-=(const Coord& op);
  //  Coord operator-=(const Coord& *op) const;
  Coord operator*=(double op);
  bool operator==(const Coord& op) const;
  bool operator!=(const Coord& op) const;
  double get_x() const {return x;};
  double get_y() const {return y;};
  double get_z() const {return z;};
  void setX(double inX){x=inX;};
  void setY(double inY){y=inY;};
  void setZ(double inZ){z=inZ;};
  void moveTo(double inX, double inY, double inZ){x=inX; y=inY; z=inZ;};
  void moveTo(Coord inC){x=inC.get_x(); y=inC.get_y(); z=inC.get_z();};
  double distTo(const Coord& inC) const;
  double norm() const;
  double innerProduct(const Coord& inC) const;
  Coord outerProduct(const Coord& inC) const;
  double angle(const Coord& a, const Coord& b, const Coord& c) const;
  double angle(const Coord& inC) const;
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
