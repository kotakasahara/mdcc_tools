#ifndef __KPRML_H__
#define __KPRML_H__

#include "Config.h"
#include "Read.h"
#include "Write.h"
#include "EMAlgorithm.h"
#include "VBgmm.h"
#include "VBfull.h"
#include "Gaussian.h"
#include <iostream>
#include <set>
#include <map>
#include <sstream>
#include <iomanip>
#include <cstdlib>
using namespace std;

class Kprml{
 private:
  Config cfg;
 public:
  Kprml();
  int setup(int argn, char* argv[]);
  int mainStream();
  int testMode();
  int emMode();
  int vbgmmMode();
  int vbfullMode();
  int calGaussianMixtureJeffreysDistanceMode();
 int calGaussianMixtureJeffreysDistance(map<string,GaussianMixture> *gm,
					ublas::vector<double> v_init,
					ublas::vector<double> v_stepwidth,
					ublas::vector<int> v_nsteps,
					bool xyzToSpherical);
 int calGaussianElementJeffreysDistanceMode();
 int calGaussianElementJeffreysDistance(map<string,GaussianMixture> *gm,
					ublas::vector<double> v_init,
					ublas::vector<double> v_stepwidth,
					ublas::vector<int> v_nsteps,
					bool xyzToSpherical);
};
#endif
