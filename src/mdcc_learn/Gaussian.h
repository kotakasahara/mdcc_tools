#ifndef __GAUSSIAN__
#define __GAUSSIAN__

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/version.hpp>
#include "Coord.h"
using namespace std;
//using namespace boost::numeric::ublas;
namespace ublas=boost::numeric::ublas;


class Gaussian{
 private:
  int id;
  string type;
  ublas::vector<double> mu;
  ublas::matrix<double> sigma;
  ublas::vector<double> mesh;
  ublas::matrix<double> sigma_inv;
  double sigma_det;
  double mesh_sum;
 public:
  Gaussian();
  Gaussian(int inId, string inType,
	   ublas::vector<double> inMu, ublas::matrix<double> inSigma);
  ublas::vector<double> getMu(){return mu;};
  ublas::matrix<double> getSigma(){return sigma;};
  double calProb(ublas::vector<double> x);
  int getId(){return id;}
  double getMesh(int i){return mesh[i];};
  double getMeshSum(){return mesh_sum;};
  double calGaussMaharanobisDist(ublas::vector<double> x);
  double calKLDivergenceMesh(Gaussian *ge);
  double calJeffreysDistanceMesh(Gaussian *ge);
  double setProbabilityMesh(ublas::vector<double> init,
			    ublas::vector<double> stepwidth,
			    ublas::vector<int> n_steps,
			    bool xyzToSpherical,
			    double laplace_delta);
  void clearProbabilityMesh(){mesh.clear();};
  double setProbabilityMeshSum(ublas::vector<double> init,
			       ublas::vector<double> stepwidth,
			       ublas::vector<int> n_steps,
			       bool xyzToSpherical,
			       double laplace_delta);
  double calKlDivergence(ublas::vector<double> init,
			 ublas::vector<double> stepwidth,
			 ublas::vector<int> n_steps,
			 bool xyzToSpherical,
			 double laplace_delta,
			 Gaussian* ge);
};

class GaussianMixture{
 private:
  string type;
  vector<double> pi;
  vector<Gaussian> gauss;
  ublas::vector<double> mesh;
  double mesh_sum;
 public:
  GaussianMixture();
  void setType(string inType){type=inType;};
  //  GaussianMixture(vector<double> inPi, vector<Gaussian> gauss);
  int pushMixtureElement(double inPi, Gaussian inElem);
  double calProb(ublas::vector<double> x);
  double calKLDivergenceMesh(GaussianMixture *gm);
  double calJeffreysDistanceMesh(GaussianMixture *gm);
  double setProbabilityMesh(ublas::vector<double> init,
			    ublas::vector<double> stepwidth,
			    ublas::vector<int> n_steps,
			    bool xyzToSpherical,
			    double laplace_delta);
  double getMesh(int i){return mesh[i];};
  double getSumMesh(){return mesh_sum;};
  vector<Gaussian>::iterator getGaussBegin(){return gauss.begin();}
  vector<Gaussian>::iterator getGaussEnd(){return gauss.end();}
};

#endif
