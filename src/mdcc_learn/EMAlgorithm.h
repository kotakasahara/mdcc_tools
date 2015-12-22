#ifndef __EM_ALGORITHM_H__
#define __EM_ALGORITHM_H__

#include "define.h"
#include "math.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/version.hpp>

//using namespace boost::numeric::ublas;
namespace ublas=boost::numeric::ublas;
//#include "math.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
//using namespace std;

class EMAlgorithm{
 private:
 protected:
  std::string fnOutGaussian;

  int nMixedElem;
  int dimension;
  //dataTable[n][i] => [data][feature]
  ublas::vector<ublas::vector<double> > dataTable;
  //dataTable[k][i] => [gaussian_element][feature]
  ublas::vector<ublas::vector<double> > mu;
  //sigma[k](i,j) => [gaussian_element](feature,feature)
  ublas::vector<ublas::matrix<double> > sigma;
  //pi[k] => [gaussian_element]
  ublas::vector<double> pi;
  //responsibility(n,k) => (data,gaussian_element)
  ublas::vector<ublas::vector<double> > responsibility;

  int maxIteration;
  int maxIterationKmeans;
  double thresholdConvergence;
  int printInterval;
  //  int maxNKmeans;
  double minDeterminantSigma;

  ublas::vector<double> sum_data;
  ublas::vector<double> mean_data;
 public:
  EMAlgorithm();
  EMAlgorithm(const ublas::vector<ublas::vector<double> > &inData,
	      std::string inFnOutGaussian,
	      int inNMixedElement, int inMaxIteration,
	      double inThresholdConvergence,
	      int inPrintInterval,
	      int inMaxIterationKmeans);
  int getDimension(){return dimension;};
  int getNMixedElem(){return nMixedElem;};
  int initialize();
  int kmeansRepeat(int maxNKmeans);
  int iteration();
  ublas::vector<ublas::vector<double> > stepE();
  int stepM(ublas::vector<ublas::vector<double> > &newMu,
	    ublas::vector<ublas::matrix<double> > &newSigma,
	    ublas::vector<double> &newPi);
  double calLogLikelihood();
  //  double calGaussian(ublas::vector<double> gX,
  //		     ublas::vector<double> gMu,
  //		     ublas::matrix<double> gSigma);
  void test();
  std::string outResultsString();
  double calMDL();
  int kmeans();
};


#endif
