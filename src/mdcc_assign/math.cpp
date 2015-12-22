#include "math.h"

const ublas::matrix<double> vector_to_matrix_transpose(const ublas::vector<double>& v){
  ublas::matrix<double> mtx(v.size(),1);
  for(int i=0;i<(int)v.size();i++) mtx(i,0)=v[i];
  return mtx;
}

const ublas::matrix<double> vector_to_matrix(const ublas::vector<double>& v){
  ublas::matrix<double> mtx(1,v.size());
  for(int i=0;i<(int)v.size();i++) mtx(0,i)=v[i];
  return mtx;
}

double cal_gaussian(const ublas::vector<double> &g_x,
		   const ublas::vector<double> &g_mu,
		   const ublas::matrix<double> &g_sigma){
  ublas::matrix<double> g_sigma_inv;
  invert(g_sigma, g_sigma_inv);
  return cal_gaussian_inv(g_x, g_mu, g_sigma_inv, determinant(g_sigma));
}

double cal_gaussian_inv(const ublas::vector<double> &g_x,
			const ublas::vector<double> &g_mu,
			const ublas::matrix<double> &g_sigma_inv,
			const double det){
  double val;
  ublas::vector<double> hensa_v=g_x-g_mu;
  ublas::matrix<double> hensa=vector_to_matrix(hensa_v);
  ublas::matrix<double> hensa_t=vector_to_matrix_transpose(hensa_v);
  //ublas::matrix<double> gSigmaInv;
  //std::cout << "calgauss 0"<<std::endl;
  //  invert(gSigma,gSigmaInv);
  //std::cout<<"calgauss 1 gSigmaInv:"<<gSigmaInv.size1()<<","<<gSigmaInv.size2()<<std::endl;
  //  std::cout<<"calgauss 2 "<<hensa<<std::endl;
  //  std::cout<<"calgauss 3 "<<hensaT<<std::endl;
  //  std::cout<<"calgauss 4 "<<gSigma<<std::endl;
  //  std::cout<<"calgauss 5 "<<gSigmaInv<<std::endl;
  ublas::matrix<double> tmp1=ublas::prod(g_sigma_inv,hensa_t);
  //    std::cout << "calgauss 6"<<std::endl;
  ublas::matrix<double> tmp6=ublas::prod(hensa,tmp1);
  //  double tmp6=calMaharanobisDistInv(gX,gMu,gSigmaInv);
  //    std::cout << "calgauss 7 "<<tmp6.size1()<<" "<<tmp6.size2()<<" "<<tmp6(0,0)<<std::endl;
  double tmp2=exp(-0.5*tmp6(0,0));
  //   std::cout << "calgauss 9 "<<tmp2<<std::endl;
  double tmp3=sqrt(det);
    //     std::cout << "calgauss 9.5 "<<determinant(gSigma)<<std::endl;
    //  std::cout << "calgauss 10 "<<tmp3<<std::endl;
  double tmp4=pow(2*PI,g_x.size()*0.5);
  //  std::cout << "calgauss 11 "<<tmp4<<std::endl;
  double tmp5=1/(tmp4*tmp3)*tmp2;
  //  std::cout << "calgauss 12 "<<tmp5<<std::endl;
  val=tmp5;
  if (!(val > 0 || val < 1.0)){
    std::cout << "prob " << val << " ";
    for(int i=0; i<g_x.size(); i++){
      std::cout << g_x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << " mu ";
    for(int i=0; i<g_mu.size(); i++){
     std::cout << g_mu[i] << " ";
    }
    std::cout << " sigma ";
    std::cout << std::endl;
    for(int i=0; i<g_x.size(); i++){
      for(int j=0; j<g_x.size(); j++){
	std::cout << g_sigma_inv(i,j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  return val;
}
double cal_maharanobis_dist_inv(const ublas::vector<double> &gX,
				const ublas::vector<double> &gMu,
				const ublas::matrix<double> &gSigmaInv){
  double val;
  ublas::vector<double> hensaV=gX-gMu;
  ublas::matrix<double> hensa=vector_to_matrix(hensaV);
  ublas::matrix<double> hensaT=vector_to_matrix_transpose(hensaV);
  ublas::matrix<double> tmp1=ublas::prod(gSigmaInv,hensaT);
  ublas::matrix<double> tmp6=ublas::prod(hensa,tmp1);
  
  //  for(int i=0;i<3;i++){
  //    for(int j=0;j<3;j++){
  //     std::cout<<  gSigmaInv(i,j) <<" ";
  //	}
  //    std::cout << std::endl;
  //  }
  return sqrt(tmp6(0,0));
}

double digamma(double z) {
  double EULER_MASCHERONI = -0.5772156649015328606065121;
  double DIGAMMA_COEF_1 = 1.0 / 12;
  double DIGAMMA_COEF_2 = 1.0 / 120;
  double DIGAMMA_COEF_3 = 1.0 / 252;
  double DIGAMMA_COEF_4 = 1.0 / 240;
  double DIGAMMA_COEF_5 = 1.0 / 132;
  double DIGAMMA_COEF_6 = 691.0 / 32760;
  double DIGAMMA_COEF_7 = 1.0 / 12;
  double DIGAMMA_LARGE = 9.5;
  double DIGAMMA_SMALL = .000001;
  
  double psi = 0;
  if (z < DIGAMMA_SMALL) {
    psi = EULER_MASCHERONI - (1 / z);
    return psi;
  }
  while (z < DIGAMMA_LARGE) {
    psi -= 1 / z;
    z++;
  }
  double invZ = 1 / z;
  double invZSquared = invZ * invZ;
  psi += log(z)
    - .5
    * invZ
    - invZSquared
    * (DIGAMMA_COEF_1 - invZSquared
       * (DIGAMMA_COEF_2 - invZSquared
	  * (DIGAMMA_COEF_3 - invZSquared
	     * (DIGAMMA_COEF_4 - invZSquared
		* (DIGAMMA_COEF_5 - invZSquared
		   * (DIGAMMA_COEF_6 - invZSquared * DIGAMMA_COEF_7))))));
  return psi;
}

const ublas::vector<double> zero_vector(int size){
  ublas::vector<double> zv(size);
  for(int i=0;i<size;i++){
    zv[i] = 0.0;
  }
  return zv;
}
const ublas::matrix<double> zero_matrix(int size1,int size2){
  ublas::matrix<double> zm(size1,size2);
  for(int i=0;i<size1;i++){
    for(int j=0;j<size2;j++){
      zm(i,j) = 0.0;
  }
  }
  return zm;
}
const ublas::matrix<double> id_matrix(int size){
  ublas::matrix<double> im(size,size);
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      im(i,j) = 0.0;
    }
    im(i,i) = 1.0;
  }
  return im;
}
const ublas::matrix<double> v2m(const ublas::vector<double>& v){
  ublas::matrix<double> mtx(v.size(),1);
  for(int i=0;i<(int)v.size();i++) mtx(i,0)=v[i];
  return mtx;
}

const ublas::matrix<double> v2mT(const ublas::vector<double>& v){
  ublas::matrix<double> mtx(1,v.size());
  for(int i=0;i<(int)v.size();i++) mtx(0,i)=v[i];
  return mtx;
}
