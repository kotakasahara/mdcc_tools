#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "define.h"
#include "Type.h"
#include <vector>
#include <set>
#include <cstdlib>
#include <iostream>
using namespace std;

class Config{
 private:
 public:
  int mode;
  string fn_cfg;
  string fn_result;
  string fn_interactions;
  string fn_gaussians;
  double assign_max_distance;
  double assign_min_prob_dens;
  int    flg_interaction_type;
  int interactions_n_begin;
  int interactions_n_end;

  vector<int> target_columns;
  vector<int> gmm_type;
  int skip_data;
  int skip_header;
  int skip_header_gaussian;

  int format_output;
  Config();
  ~Config();

  void set_all(int argn,char* argv[]);
  void set_all(const vector<string>& arg);
  void operator=(Config op);
};

#endif
