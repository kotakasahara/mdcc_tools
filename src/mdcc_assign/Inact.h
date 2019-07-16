#ifndef __INACT_H___
#define __INACT_H___

#include <string>
#include <sstream>
#include "Coord.h"
#include "Type.h"
using namespace std;

class Inact {
 private:
  string code;
  int id;
  vector<int> type;
 public:
  Inact(string in_code, int in_id,
	const vector<int>& in_type);
  string get_code() const {return code;};
  int get_id() const {return id;};
  const vector<int>& get_type() const {return type;};
};

class TriTriInact : public Inact{
 private:
  int ptriid;
  int residv;
  int restype;
  int ltriid;
  Coord ig1;
  Coord ig2;
  Coord ig3;
  //  Coord ig_cent;
  //  Coord ig_norm;
 public:
  TriTriInact(string in_code, int in_id,
	      const vector<int>& in_type,
	      int in_ptriid,
	      int in_resid, int in_restype,
	      int in_ltriid,
	      Coord in_ig1, Coord in_ig2, Coord in_ig3);
  int get_ptri_id() const {return ptriid;};
  int get_res_idv() const {return residv;};
  int get_ltri_id() const {return ltriid;};
  Coord get_ig1() const {return ig1;};
  Coord get_ig2() const {return ig2;};
  Coord get_ig3() const {return ig3;};
  //  const Coord& get_ig_cent() const {return ig_cent;};
  //  const Coord& get_ig_norm() const {return ig_norm;};
  //  void set_ig_cent_norm();
};

class TriAtomInact : public Inact{
 private:
  int ptriid;
  int residv;
  int restype;
  int latomidv;
  Coord ig;
  //  Coord ig_cent;
  //  Coord ig_norm;
 public:
 TriAtomInact(string in_code, int in_id,
	      const vector<int>& in_type,
	      int in_ptriid,
	      int in_resid, int in_restype,
	      int in_latomidv,
	      Coord in_ig);
  int get_ptri_id() const {return ptriid;};
  int get_res_idv() const {return residv;};
  int get_latom_idv() const {return latomidv;};
  Coord get_ig() const {return ig;};
  //  const Coord& get_ig_cent() const {return ig_cent;};
  //  const Coord& get_ig_norm() const {return ig_norm;};
  //  void set_ig_cent_norm();
};

class ResAtomInact : public Inact{
 private:
  int residv;
  int restype;
  int latomidv;
  Coord ig;
  //  Coord ig_cent;
  //  Coord ig_norm;
 public:
  ResAtomInact(string in_code, int in_id,
	      const vector<int>& in_type,
	      int in_resid, int in_restype,
	      int in_latomidv,
	      Coord in_ig);
  int get_res_idv() const {return residv;};
  int get_latom_idv() const {return latomidv;};
  Coord get_ig() const {return ig;};
  //  const Coord& get_ig_cent() const {return ig_cent;};
  //  const Coord& get_ig_norm() const {return ig_norm;};
  //  void set_ig_cent_norm();
};


#endif
