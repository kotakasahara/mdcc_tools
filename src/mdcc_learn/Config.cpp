#include "define.h"
#include "Config.h"
#include <iostream>
using namespace std;

Config::Config(){
  mode=M_VBFULL;
  fnCfg="";
  fnGaussianMixtureInteractions = "";
  klCalcPairBegin = -1;
  klCalcPairEnd = -1;
  klCalcPairBegin1 = -1;
  klCalcPairEnd1 = -1;
  klCalcPairBegin2 = -1;
  klCalcPairEnd2 = -1;
  klCalcLaplaceDelta = 0.00001;
  flgKlGridConvSpherical = false;
  maxDistanceForGaussianMixtures = 0.0;
  modeSaveMemory = 0;
  flgSameGm = false;
  dataSkip=1;
  headerSkip=0;
  endFrame=-1;
  formatDataTable = DATA_TABLE;

  maxIteration=10000;
  maxIterationKmeans=100;
  thresholdConvergence=0.0001;
  printInterval=0;
  vbfull_alpha_0 = 0.001;
  vbfull_w_0_coef = 1.0;
  vbfull_ny_0 = 5.0;
  formatDataTable=DATA_KKTRJ;

}

Config::~Config(){
}

void Config::setAll(int argn, char* argv[]){
  vector<string> arg;
  int i;
  for(i=1;i<argn;i++)
    arg.push_back(string(argv[i]));
  setAll(arg);
}
void Config::setAll(vector<string> arg){
  vector<string>::iterator itr;
  string type,val;
  for(itr=arg.begin(); itr!=arg.end(); itr++){
    if(*itr=="-mode"){
      itr++;
      if(*itr=="test")            { mode=M_TEST; }
      else if(*itr=="em")    { mode=M_EM; }
      else if(*itr=="vbgmm")    { mode=M_VBGMM; }
      else if(*itr=="vbfull")    { mode=M_VBFULL; }
      else if(*itr=="cal-gm-dist")    { mode=M_CAL_GM_DIST; }
      else if(*itr=="cal-ge-dist")    { mode=M_CAL_GE_DIST; }
      else{
	mode=M_VBFULL;
	//cerr<<"invalid mode ["<<(*itr)<<"]\n"; exit(1);
      }
    }
    else if(*itr=="-fn-cfg"){ fnCfg=*++itr; }
    else if(*itr=="-fn-data-table"){ fnDataTable=*++itr; }
    else if(*itr=="-fn-out-gaussian"){ fnOutGaussian=*++itr; }
    else if(*itr=="-format-data-table"){
      itr++;
      if(*itr=="kktrajtrans") { formatDataTable=DATA_KKTRJ; }
      else{ formatDataTable=DATA_TABLE; }
    }
    else if(*itr=="-feature"){
      int type=atoi((*++itr).c_str());
      int column=atoi((*++itr).c_str());
      string header=(*++itr);
      featureDef.push_back(FeatureDefinition(type,column,header));
    }
    else if(*itr=="-atom"){
      int type=1;
      int column=atoi((*++itr).c_str());
      string header("x");
      featureDef.push_back(FeatureDefinition(type,column,header));
    }
    //em algorithm
    else if(*itr=="-n-mixed-element"){ nMixedElement=atoi((*++itr).c_str()); }
    else if(*itr=="-n-max-iteration"){ maxIteration=atoi((*++itr).c_str()); }
    else if(*itr=="-n-max-iteration-kmeans"){ maxIterationKmeans=atoi((*++itr).c_str()); }
    else if(*itr=="-threshold-convergence"){ thresholdConvergence=atof((*++itr).c_str()); }
    else if(*itr=="-print-interval"){ printInterval=atoi((*++itr).c_str()); }
    //vbgmm, vbfull
    else if(*itr=="-vbgmm-hp-beta"){ vbgmm_hp_beta=atof((*++itr).c_str()); }
    else if(*itr=="-vbfull-alpha-0"){ vbfull_alpha_0=atof((*++itr).c_str()); }
    else if(*itr=="-vbfull-w-0-coef"){ vbfull_w_0_coef=atof((*++itr).c_str()); }
    else if(*itr=="-vbfull-ny-0"){ vbfull_ny_0=atof((*++itr).c_str()); }
    // cal Gaussian Mixture Jeffreys Distance Mode
    else if(*itr=="-fn-gaussianmixture"){ fnGaussianMixtureInteractions=*++itr; }
    else if(*itr=="-gaussianmixture-type-length"){ gaussianMixtureTypeLength=atoi((*++itr).c_str()); }
    else if(*itr=="-max-distance-for-gaussianmixtures"){ maxDistanceForGaussianMixtures=atof((*++itr).c_str()); }
    else if(*itr=="-kl-grid-init"){ klGridInit.push_back(atof((*++itr).c_str())); }
    else if(*itr=="-kl-grid-stepwidth"){ klGridStepwidth.push_back(atof((*++itr).c_str())); }
    else if(*itr=="-kl-grid-n-steps"){ klGridNSteps.push_back(atoi((*++itr).c_str())); }
    else if(*itr=="-kl-calc-pair-begin"){ klCalcPairBegin=atoi((*++itr).c_str()); }
    else if(*itr=="-kl-calc-pair-end"){ klCalcPairEnd=atoi((*++itr).c_str()); }
    else if(*itr=="-kl-calc-pair-begin-1"){ klCalcPairBegin1=atoi((*++itr).c_str()); }
    else if(*itr=="-kl-calc-pair-end-1"){ klCalcPairEnd1=atoi((*++itr).c_str()); }
    else if(*itr=="-kl-calc-pair-begin-2"){ klCalcPairBegin2=atoi((*++itr).c_str()); }
    else if(*itr=="-kl-calc-pair-end-2"){ klCalcPairEnd2=atoi((*++itr).c_str()); }
    else if(*itr=="-kl-calc-laplace-delta"){ klCalcLaplaceDelta=atof((*++itr).c_str()); }
    else if(*itr=="-flg-kl-grid-conv-spherical"){ flgKlGridConvSpherical=true; }
    else if(*itr=="-kl-max-maharanobis-distance"){ klMaxMaharanobisDist=atof((*++itr).c_str()); }
    else if(*itr=="-mode-save-memory"){ modeSaveMemory=atoi((*++itr).c_str()); }
    else if(*itr=="-skip-data"){ dataSkip=atoi((*++itr).c_str()); }
    else if(*itr=="-skip-header"){ headerSkip=atoi((*++itr).c_str()); }
    else if(*itr=="-end-frame"){ endFrame=atoi((*++itr).c_str()); }
    else if(*itr=="-flg-same-gm"){ flgSameGm=true; }
    else{
      cerr<<"unknown keyword <"<<(*itr)<<">"<<endl;
    }
  }
}

