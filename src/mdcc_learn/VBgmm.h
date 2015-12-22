#ifndef __VB_GMM_H__
#define __VB_GMM_H__
#include "EMAlgorithm.h"
class VBgmm : public EMAlgorithm {
 private:

  //Hyper-parameters
  double hp_beta;
  double hp_mean_g;
  ublas::matrix<double> hp_prec_g;
  double hp_dof_w;
  ublas::matrix<double> hp_scale_w;

  //prior
  //p (->s) //[component,data] 
  ublas::vector<ublas::vector<double> > prior_s;
  //m (->mu) Gauss //[component,dimension(*dimension)]
  ublas::vector<ublas::vector<double> > mean_prior_mean;        //m
  ublas::vector<ublas::matrix<double> > prec_prior_mean;     //T
  //sigma-Wishart  //[component,dimension(*dimension)]
  ublas::vector<double> dof_prior_prec;     //ny
  ublas::vector<ublas::matrix<double> > scale_prior_prec;  //V

 public:
  VBgmm();
  VBgmm(const ublas::vector<ublas::vector<double> > &inData,
	std::string inFnOutGaussian,
	int inNMixedElement,
	int inMaxIteration,
	double inThresholdConvergence,
	int inPrintInterval,
	double in_hp_beta,
	int inMaxIterationKmeans);
  int initialize();
  int iteration();
  int stepE();
  int stepE_prior_for_mean();
  int stepE_prior_for_prec();
  int stepM();
  double cal_lower_bound();
  int main();
  int set_ml_parameters();
  std::string get_results_st();
};

#endif
