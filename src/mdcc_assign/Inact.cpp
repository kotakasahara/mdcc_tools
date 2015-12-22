#include "Inact.h"

Inact::Inact(string in_code,
	     int in_id,
	     const vector<int>& in_type){
  code = in_code;
  id = in_id;
  type = in_type;
}

TriTriInact::TriTriInact(string in_code, int in_id,
			 const vector<int>& in_type,
			 int in_ptriid,
			 int in_resid,
			 int in_restype,
			 int in_ltriid,
			 Coord in_ig1, Coord in_ig2, Coord in_ig3) :
  Inact(in_code, in_id, in_type){
  
  ptriid = in_ptriid;
  residv = in_resid;
  restype = in_restype;
  ltriid = in_ltriid;
  ig1 = in_ig1;
  ig2 = in_ig2;
  ig3 = in_ig3;
}
TriAtomInact::TriAtomInact(string in_code, int in_id,
			   const vector<int>& in_type,
			   int in_ptriid,
			   int in_resid,
			   int in_restype,
			   int in_latomidv,
			   Coord in_ig) :
  Inact(in_code, in_id, in_type){
  ptriid = in_ptriid;
  residv = in_resid;
  restype = in_restype;
  latomidv = in_latomidv;
  ig = in_ig;
}
