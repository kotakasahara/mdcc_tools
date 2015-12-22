#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "define.h"
#include "FeatureDefinition.h"
#include <vector>
#include <set>
#include <cstdlib>
using namespace std;

class Config{
 private:
 public:
  int mode;
  string fnCfg;
  string fnDataTable;
  string fnOutGaussian;
  int formatDataTable;
  int nMixedElement;
  int maxIteration;
  int maxIterationKmeans;
  double thresholdConvergence;  
  int printInterval;
  double vbgmm_hp_beta;
  double vbfull_alpha_0;
  double vbfull_w_0_coef;
  double vbfull_ny_0;
  vector<FeatureDefinition> featureDef;
  string fnGaussianMixtureInteractions;
  int gaussianMixtureTypeLength;
  double maxDistanceForGaussianMixtures;
  vector<double> klGridStepwidth;
  vector<double> klGridInit;
  vector<int> klGridNSteps;
  int klCalcPairBegin;
  int klCalcPairEnd;
  int klCalcPairBegin1;
  int klCalcPairEnd1;
  int klCalcPairBegin2;
  int klCalcPairEnd2;
  double klCalcLaplaceDelta;
  bool flgKlGridConvSpherical;
  bool flgSameGm;
  int modeSaveMemory;
  int dataSkip;
  int headerSkip;
  int endFrame;
  double klMaxMaharanobisDist;
  Config();
  ~Config();

  void setAll(int argn,char* argv[]);
  void setAll(vector<string> arg);
  void operator=(Config op);
};

#endif
