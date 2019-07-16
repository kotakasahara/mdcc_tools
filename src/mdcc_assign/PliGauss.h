#ifndef __PLI_GAUSS_H__
#define __PLI_GAUSS_H__

#include "define.h"
#include "Config.h"
#include "Read.h"
#include "Write.h"

#include <iostream>
using namespace std;

class PliGauss{
 private:
  Config cfg;
  int n_row;
 public:
  PliGauss();
  ~PliGauss();
  int setup(int argn, char* argv[]);
  int mainStream();
  int test();
  int assign_tti();
  int assign_tai();
  int assign_rai();
  int assign_datatable();
  int assign_trajtrans();
};

#endif
