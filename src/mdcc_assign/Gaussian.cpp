#include "Gaussian.h"
#include "math.h"

/////////////////// Gaussian /////////////////////

Gaussian::Gaussian(){
}
Gaussian::Gaussian(int in_id, const vector<int>& in_type,
		   const ublas::vector<double>& in_mu,
		   const ublas::matrix<double>& in_sigma){
  id = in_id;
  type = in_type;
  mu = in_mu;
  sigma = in_sigma;
  invert(sigma, sigma_inv);
  sigma_det = determinant(sigma);
}
double Gaussian::cal_prob(const ublas::vector<double>& x) const {
  return cal_gaussian_inv(x, mu, sigma_inv, sigma_det);
}
double Gaussian::cal_prob(const vector<double>& x) const {
  ublas::vector<double> x_ublas(x.size());
  for(int i = 0; i < x.size(); i++){
    x_ublas[i] = x[i];
  }
  return cal_gaussian_inv(x_ublas, mu, sigma_inv, sigma_det);
}

double Gaussian::cal_gauss_maharanobis_dist(const ublas::vector<double>& x) const {
  //cout << (int)x.size()  << " " << (int)mu.size() <<" " << (int)sigma_inv.size1() << <<" " <<(int)sigma_inv.size2() << endl;;
  return cal_maharanobis_dist_inv(x, mu, sigma_inv);
}
double Gaussian::cal_gauss_maharanobis_dist(const vector<double>& x) const {
  ublas::vector<double> x_ublas(x.size());
  for(int i = 0; i < x.size(); i++){
    x_ublas[i] = x[i];
  }
  return cal_maharanobis_dist_inv(x_ublas, mu, sigma_inv);
}

/////////////////// GaussianMixture /////////////////////

GaussianMixture::GaussianMixture(){
}
int GaussianMixture::push_mixture_element(double in_pi, const Gaussian& in_elem){
  pi.push_back(in_pi);
  gauss.push_back(in_elem);
  return gauss.size();
}

double GaussianMixture::cal_prob(const ublas::vector<double>& x) const {
  vector<double>::const_iterator itr_pi = pi.begin();
  vector<Gaussian>::const_iterator itr_g = gauss.begin();
  double p = 0;
  for(;itr_pi!=pi.end();itr_pi++, itr_g++){
    p += (*itr_pi)*(*itr_g).cal_prob(x);
  }
  return p;
}

