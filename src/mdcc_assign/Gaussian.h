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
  vector<int> type;
  ublas::vector<double> mu;
  ublas::matrix<double> sigma;
  ublas::vector<double> mesh;
  ublas::matrix<double> sigma_inv;
  double sigma_det;
  double mesh_sum;
 public:
  Gaussian();
  Gaussian(int in_id, const vector<int>& in_type,
	   const ublas::vector<double>& in_mu,
	   const ublas::matrix<double>& in_sigma);
  ublas::vector<double> get_mu() const {return mu;};
  ublas::matrix<double> get_sigma() const {return sigma;};
  double cal_prob(const ublas::vector<double>& x) const;
  double cal_prob(const vector<double>& x) const;
  int get_id() const {return id;}
  double get_mesh(int i) const {return mesh[i];};
  double get_mesh_sum() const {return mesh_sum;};
  double cal_gauss_maharanobis_dist(const ublas::vector<double>& x) const;
  double cal_gauss_maharanobis_dist(const vector<double>& x) const;
};

class GaussianMixture{
 private:
  vector<int> type;
  vector<double> pi;
  vector<Gaussian> gauss;
  ublas::vector<double> mesh;
  double mesh_sum;
 public:
  GaussianMixture();
  void set_type(const vector<int>& in_type) {type=in_type;};
  //  GaussianMixture(vector<double> inPi, vector<Gaussian> gauss);
  int push_mixture_element(double in_pi, const Gaussian& in_elem);
  double cal_prob(const ublas::vector<double>& x) const;
  const vector<Gaussian>& get_gauss() const {return gauss;};
  const vector<double>& get_pi() const {return pi;};
};

#endif
