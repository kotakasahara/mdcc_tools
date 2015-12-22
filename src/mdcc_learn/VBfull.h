#ifndef __VB_FULL_H__
#define __VB_FULL_H__

#include <cmath>
#include "EMAlgorithm.h"
#include "Write.h"
using namespace std;

class VBfull : public EMAlgorithm {
 private:
  //prior
  //Dirichlet
  double alpha_0; 
  //Gauss
  //[dimension]
  ublas::vector<double> m_0;
  double beta_0;
  //Wishart
  //[dimension, dimension]
  ublas::matrix<double> W_0;
  double ny_0;
  
  //
  // [mixed elem][dimension]
  ublas::vector<double> alpha_k;
  // [mixed elem][dimension]
  ublas::vector<ublas::vector<double> > m_k;
  // [mixed elem]
  ublas::vector<double> beta_k;
  // [mixed elem][dimension,dimension]
  ublas::vector<ublas::matrix<double> > W_k;
  // [mixed elem]
  ublas::vector<double> ny_k;

  // [mixed elem]
  ublas::vector<double> N_k;
  // [mixed elem][dimension]
  ublas::vector<ublas::vector<double> > x_bar_k;
  // [mixed elem][dimension,dimension]
  ublas::vector<ublas::matrix<double> > S_k;

  // 
  ublas::vector<double> ln_lambda_til_k;
  ublas::vector<double> ln_pi_til_k;

  //[data][mixed elem]
  //  ublas::vector<ublas::vector<double> > rho;

 public:
  VBfull();
  VBfull(const ublas::vector<ublas::vector<double> > &inData,
	 string inFnOutGaussian,
	 int inNMixedElement,
	 int inMaxIteration,
	 double inThresholdConvergence,
	 int inPrintInterval,
	 double in_alpha_0,
	 double in_w_0_coef,
	 double in_ny_0,
	 int inMaxIterationKmeans);
  int initialize(double in_alpha_0,
		 double in_w_0_coef,
		 double in_ny_0);
  int update_N_xbar_S();
  int update_alpha_beta_ny_m_W();
  int update_ln_lambda_pi();
  int update_responsibility();
  double cal_lower_bound();
  int main();
  int iteration();
  int set_ml_parameters();
  std::string get_results_st();
  std::string get_results_st_tsv();
};
#endif
