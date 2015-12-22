#ifndef __WRITE_H__
#define __WRITE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "Inact.h"
using namespace std;

class Write{
 private:
  string filename;
  bool op;
 public:
  ofstream ofs;
  Write(string inFn);
  string getFn(){return filename;};
  bool is_open(){return op;};
  int open();
  int openApp();
  int close();
  int write_assigned_tti();
  int write_line(const string& st);
  int write_assign_data(double** assign, int n_gauss, int n_row);
};

#endif
