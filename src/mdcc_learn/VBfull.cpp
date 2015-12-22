#include "VBfull.h"
#include "math.h"
VBfull::VBfull() : EMAlgorithm(){
}

VBfull::VBfull(const ublas::vector<ublas::vector<double> > &inData,
	       std::string inFnOutGaussian,
	       int inNMixedElement,
	       int inMaxIteration,
	       double inThresholdConvergence,
	       int inPrintInterval,
	       double in_alpha_0,
	       double in_w_0_coef,
	       double in_ny_0,
	       int inMaxIterationKmeans)
  : EMAlgorithm(inData, inFnOutGaussian, 
		inNMixedElement, inMaxIteration,
		inThresholdConvergence, inPrintInterval,
		inMaxIterationKmeans){
  initialize(in_alpha_0, in_w_0_coef, in_ny_0);
}

int VBfull::initialize(double in_alpha_0, double in_w_0_coef, double in_ny_0){
  EMAlgorithm::initialize();
  
  alpha_0 = in_alpha_0;
  m_0 = zero_vector(dimension);
  beta_0 = 1.0;
  W_0 = id_matrix(dimension) * in_w_0_coef;
  ny_0 = in_ny_0;

  alpha_k = zero_vector(nMixedElem);
  m_k = ublas::vector<ublas::vector<double> >(nMixedElem);
  beta_k = zero_vector(nMixedElem);
  W_k = ublas::vector<ublas::matrix<double> >(nMixedElem);
  ny_k = zero_vector(nMixedElem);
  N_k = zero_vector(nMixedElem);
  x_bar_k = ublas::vector<ublas::vector<double> >(nMixedElem);
  S_k = ublas::vector<ublas::matrix<double> >(nMixedElem);
  for(int k=0; k<nMixedElem; k++){
    m_k[k] = zero_vector(dimension);
    W_k[k] = W_0;
    x_bar_k[k] = zero_vector(dimension);
    S_k[k] = W_0;
    ny_k[k] = ny_0;
    alpha_k[k] = alpha_0;
    beta_k[k] = beta_0;
  }
  ln_pi_til_k = zero_vector(nMixedElem);
  ln_lambda_til_k = zero_vector(nMixedElem);
  sum_data = ublas::vector<double>(dimension);
  mean_data = ublas::vector<double>(dimension);
  sum_data = zero_vector(dimension);
  mean_data = zero_vector(dimension);
  for(int n=0; n<(int)dataTable.size(); n++){
    for(int d=0; d<dimension; d++){
      sum_data[d] += dataTable[n][d];
    }
  }

  for(int d=0; d<dimension; d++){
    mean_data[d] = sum_data[d] / (double)dataTable.size();
  }
  for(int k=0; k<nMixedElem; k++){
    m_k[k] = mean_data;
  }
  for(int n=0; n<(int)dataTable.size(); n++){
    for(int k=0; k<nMixedElem; k++){
      responsibility[n][k] = 1.0/(dataTable.size()*nMixedElem);
    }
  }
  //  cout << "init n_xbar_s"<<endl;
  //  update_N_xbar_S();
  //  cout << "init ln_lambda_pi"<<endl;
  //  update_ln_lambda_pi();
  return 0;
}

int VBfull::update_N_xbar_S(){
  for(int k=0;k<nMixedElem;k++){
    N_k[k] = 0.0;
    x_bar_k[k] = zero_vector(dimension);
    for(int n=0;n<(int)dataTable.size();n++){
      N_k[k] += responsibility[n][k];
      x_bar_k[k] += responsibility[n][k] * dataTable[n];
    }
    if(N_k[k] > 0){
      x_bar_k[k] /= N_k[k];
    }else{
      x_bar_k[k] = zero_vector(dimension);
    }
    //cout << " k:"<<k<<" n_k[k]="<<N_k[k]<<" x_bar_k[k]="<<x_bar_k[k]<<endl;
    S_k[k] = zero_matrix(dimension,dimension);
    if(N_k[k] > 0){
      for(int n=0;n<(int)dataTable.size();n++){
	ublas::vector<double> v_xxb = dataTable[n]-x_bar_k[k];
	S_k[k] += responsibility[n][k]*ublas::prod(v2m(v_xxb),v2mT(v_xxb));
      }
      S_k[k] /= N_k[k];
    }else{
      S_k[k] = W_0;
    }
    //    cout << "N_k["<<k<<"] "<<N_k[k]<<endl;
  }
  return 0;
}

int VBfull::update_alpha_beta_ny_m_W(){
  for(int k=0;k<nMixedElem;k++){
    //    cout << "dbg2 "<<k<<endl;
    alpha_k[k] = alpha_0 + N_k[k];
    //    cout << "dbg2a "<<alpha_k[k]<<" "<<alpha_0<<" "<<N_k[k]<<endl;
    beta_k[k] = beta_0 + N_k[k];
    //    cout << "dbg2b "<<beta_k[k]<<endl;
    ny_k[k] = ny_0 + N_k[k];
    //    cout << "dbg2n "<<ny_k[k]<<endl;
    m_k[k] = 1/beta_k[k] * (beta_0*m_0 + N_k[k]*x_bar_k[k]);
    //cout << "dbg2m "<<m_k[k][0]<<endl;
    ublas::vector<double> v_xbm = x_bar_k[k] - m_0;
    //    cout << "x_bar_k[k] "<<x_bar_k[k]<<endl;
    ublas::matrix<double> W_0_inv;
    //    cout << "dbg2 1"<<endl;
    invert(W_0,W_0_inv);
    //    cout << "dbg2 2"<<endl;
    invert(W_0_inv + N_k[k]*S_k[k]
	   + (beta_0*N_k[k])/(beta_0+N_k[k])
	   * (ublas::prod(v2m(v_xbm),v2mT(v_xbm))),
	   W_k[k]);
    //    cout << "dbg2 3"<<endl;
    //    cout << "N_k[k] : "<<N_k[k] << endl;
    //    cout << "S_k[k] : "<<S_k[k] << endl;
    //    cout << "prod(v_xbm,v_xbmT) "<< ublas::prod(v2m(v_xbm),v2mT(v_xbm)) << endl;

  }
  return 0;
}

int VBfull::update_ln_lambda_pi(){
  double alpha_hat=0;
  for(int k=0;k<nMixedElem;k++){
    alpha_hat += alpha_k[k];
  }
  double digamma_alpha_hat = digamma(alpha_hat);
  double dim_log2 = dimension * log(2);
  for(int k=0;k<nMixedElem;k++){
    ln_lambda_til_k[k] = 0;
    for(int i=1;i<=dimension;i++){
      ln_lambda_til_k[k] += digamma((ny_k[k] + 1 - i)*0.5);
    }
    //    cout << "ln_lambda_til_k["<<k<<"] 1 "<<ln_lambda_til_k[k]<<endl;
    //    cout << "det W_k " << determinant(W_k[k]) << endl;
    //    cout << "W_k[k] "<<W_k[k]<<endl;
    ln_lambda_til_k[k] += dim_log2 + log(determinant(W_k[k]));
    ln_pi_til_k[k] = digamma(alpha_k[k]) - digamma_alpha_hat;
    //    cout << "ln_lambda_til_k["<<k<<"] "<<ln_lambda_til_k[k]<<endl;
    //    cout << "ln_pi_til_k["<<k<<"] "<<ln_pi_til_k[k]<<endl;
    //        cout << "digamma(alpha_k[k]) "<<digamma(alpha_k[k])
    //    	 <<" alpha_k[k] " << alpha_k[k]
    //    	 <<" digamma(alpha_hat) " <<digamma_alpha_hat
    //    	 <<" alpha_hat "<<alpha_hat<<endl;
  }
  return 0;
}

int VBfull::update_responsibility(){
  //  cout << "dbg 1"<<endl;
  ublas::vector<ublas::vector<double> > rho_n_k(dataTable.size());
  //  cout << "dbg 2"<<endl;
  for(int n=0;n<(int)dataTable.size();n++){  
    rho_n_k[n] = ublas::vector<double>(nMixedElem);
  }
  //  cout << "dbg 3"<<endl;
  double d_ln_2pi = dimension*0.5 * log(PI+PI);
  for(int k=0;k<nMixedElem;k++){
    //        cout << "dbg 4 "<<k<<endl;
    for(int n=0;n<(int)dataTable.size();n++){  
      ublas::vector<double> xm;
      //      cout << "dbg 4.1 "<<endl;
      xm = dataTable[n] - m_k[k];
      //      cout << "dbg 4.2 "<<xm<<endl;
      ublas::matrix<double> tmp_prod1 = ublas::prod(v2mT(xm),W_k[k]);
      //      cout << "dbg 4.3 "<<tmp_prod1<<endl;
      double tmp_prod2 = ublas::prod(tmp_prod1,v2m(xm))(0,0);
      //      cout << "dbg 4.4 "<<tmp_prod2<<endl;
      double ex_xmu_lambda_xmu = dimension / beta_k[k] + ny_k[k] * tmp_prod2;
      //      cout << "dbg 4.5 "<<ex_xmu_lambda_xmu<<endl;
      //      cout << "dbg 4.6 "<<ln_pi_til_k[k]<<endl;
      //      cout << "dbg 4.7 "<<ln_lambda_til_k[k]<<endl;
      rho_n_k[n][k] = exp(ln_pi_til_k[k] + 0.5*ln_lambda_til_k[k] - d_ln_2pi
			  -0.5 * ex_xmu_lambda_xmu);
      //      cout << "dbg 4.8 "<<ln_pi_til_k[k] + 0.5*ln_lambda_til_k[k] - d_ln_2pi -0.5 * ex_xmu_lambda_xmu<<endl;
      //      cout << "dbg 4.9 "<<n<<" "<<k<<" "<<rho_n_k[n][k]<<endl;
    }
    //    cout << "dbg 5"<<k<<endl;
  }
  //  cout << "dbg 6"<<endl;
  ublas::vector<double> rho_n(dataTable.size());
  for(int n=0;n<(int)dataTable.size();n++){
    rho_n[n] = 0.0;
    for(int k=0;k<nMixedElem;k++){  
      rho_n[n] += rho_n_k[n][k];
    }
  }
  //  cout << "dbg 7"<<endl;
  //  double tmp=0;
  for(int k=0;k<nMixedElem;k++){
    for(int n=0;n<(int)dataTable.size();n++){  
      //      cout << "dbg 8 "<<n<<" "<<k<<" "<<rho_n_k[n][k]<<" "<<rho_n[n]<<endl;
      responsibility[n][k] = rho_n_k[n][k] / rho_n[n];
      //      tmp += responsibility[n][k];
      //                  cout << "dbg 8.1 "<<responsibility[n][k]<<endl;
    }
  } 
  
  //    cout << "dbg 9 "<<tmp<<endl;
  return 0;
}

double VBfull::cal_lower_bound(){
  //eq 10.71
  double tmp1=0;
  double dln2pi = -dimension*log(PI+PI);
  for(int k=0;k<nMixedElem;k++){
    ublas::vector<double> xbm = x_bar_k[k]-m_k[k];
    ublas::matrix<double> tmp_prod1 = ublas::prod(v2mT(xbm),W_k[k]);
    double tmp_prod2 = ublas::prod(tmp_prod1, v2m(xbm))(0,0);
    tmp1 += N_k[k] *
      ( ln_lambda_til_k[k] - dimension/beta_k[k]
	- ny_k[k] * trace(ublas::prod(S_k[k],W_k[k]))
	- ny_k[k] * tmp_prod2
	+ dln2pi
	);
  }
  double ex_ln_p_x = 0.5*tmp1;
  //eq 10.72
  double ex_ln_p_z = 0;
  for(int n=0;n<(int)dataTable.size();n++){  
    for(int k=0;k<nMixedElem;k++){
      ex_ln_p_z += responsibility[n][k]*ln_pi_til_k[k];
      //      cout <<"n:"<<n<<" k:"<<k<<" resp:"<<responsibility[n][k]
      //	   <<" ln_pi_til_k:"<<ln_pi_til_k[k]<<endl;
    }
  }
  //eq 10.73
  double ex_ln_p_pi=0;
  double sum_lgamma_alpha_0 = nMixedElem * lgamma(alpha_0);
  double lgamma_sum_alpha_0 = lgamma(alpha_0*nMixedElem);
  double ln_C_alpha = lgamma_sum_alpha_0 - sum_lgamma_alpha_0;
  double sum_ln_pi_til = 0;
  for(int k=0;k<nMixedElem;k++){
    sum_ln_pi_til += ln_pi_til_k[k];
  }
  ex_ln_p_pi = ln_C_alpha+(alpha_0-1)*sum_ln_pi_til;
  //eq 10.74
  double ex_ln_p_mu_lambda = 0;
  double tmp2=0;
  double d_ln_beta_2pi = dimension * log(beta_0/(PI+PI));
  for(int k=0;k<nMixedElem;k++){
    ublas::vector<double> mkm0 = m_k[k]-m_0;
    ublas::matrix<double> tmp_prod1 = ublas::prod(v2mT(mkm0),W_k[k]);
    double tmp_prod2 = ublas::prod(tmp_prod1,v2m(mkm0))(0,0);
    tmp2 += d_ln_beta_2pi + ln_lambda_til_k[k]
      - dimension*beta_0/beta_k[k]
      - beta_0*ny_k[k] * tmp_prod2;
  }
  tmp2 *= 0.5;
  double tmp3;
  double tmp_lgam=0.0;
  for(int i=1;i<=dimension;i++){
    tmp_lgam += lgamma((ny_0+1-i)*0.5);
  }
  tmp3 = -ny_0*0.5*log(determinant(W_0))
    - ( ny_0*dimension*0.5 * log(2)
    + dimension*(dimension-1)*0.25 * log(PI)
	+ tmp_lgam);
  double tmp4=0;
  for(int k=0;k<nMixedElem;k++){
    tmp4 += ln_lambda_til_k[k];
  }
  tmp4 *= (ny_0-dimension-1)*0.5;
  double tmp5=0;
  for(int k=0;k<nMixedElem;k++){
    ublas::matrix<double> W_0_inv;
    invert(W_0,W_0_inv);
    tmp5 += ny_k[k]*trace(ublas::prod(W_0_inv,W_k[k]));
  }
  tmp5 *= -0.5;
  ex_ln_p_mu_lambda = tmp2+nMixedElem*tmp3+tmp4+tmp5;
  //eq 10.75
  double ex_q_z = 0;
  for(int n=0;n<(int)dataTable.size();n++){  
    for(int k=0;k<nMixedElem;k++){
      if(responsibility[n][k] > 0){
	ex_q_z += responsibility[n][k]  * log(responsibility[n][k]);
      }
    }
  }  
  //eq 10.76
  double ex_q_pi = 0;
  double sum_alpha = 0;
  double sum_lgamma_alpha = 0;
  for(int k=0;k<nMixedElem;k++){
    sum_alpha += alpha_k[k];
    sum_lgamma_alpha += lgamma(alpha_k[k]);
  }
  double lgamma_alpha_hat = lgamma(sum_alpha);
  double C_alpha_k = lgamma_alpha_hat - sum_lgamma_alpha;
  double tmp6 = 0;
  for(int k=0;k<nMixedElem;k++){
    tmp6 += (alpha_k[k]-1)*ln_pi_til_k[k];
    //    cout << "k:"<<k<<" "<<tmp6<<endl;
  }
  ex_q_pi = tmp6 + C_alpha_k;
  //eq 10.77
  double ex_q_mu_lambda = 0;
  double dim2 = dimension*0.5;
  for(int k=0;k<nMixedElem;k++){
    double tmp_H = 0;
    double tmp_ln_B;
    tmp_lgam=0.0;
    for(int i=1;i<=dimension;i++){
      tmp_lgam += lgamma((ny_k[k]+1-i)*0.5);
    }
    tmp_ln_B = -ny_0*0.5*log(determinant(W_k[k]))
      - ( ny_k[k]*dimension*0.5 * log(2)
	  + dimension*(dimension-1)*0.25 * log(PI)
	  + tmp_lgam);
    tmp_H = - tmp_ln_B - (ny_k[k]-dimension-1)*0.5
      * ln_lambda_til_k[k]
      + (ny_k[k]*dimension)*0.5;
    ex_q_mu_lambda += 0.5 * ln_lambda_til_k[k]
      + dim2 * log(beta_k[k]/(PI+PI)) - dim2 - tmp_H;
  }
  cout<< "lower "<<ex_ln_p_x<<" "
      << ex_ln_p_z <<" "
      << ex_ln_p_pi << " "
      << ex_ln_p_mu_lambda <<" "
      << - ex_q_z <<" "
      << - ex_q_pi<<" "
      << - ex_q_mu_lambda<<endl;
  return ex_ln_p_x + ex_ln_p_z + ex_ln_p_pi + ex_ln_p_mu_lambda
    - ex_q_z - ex_q_pi - ex_q_mu_lambda;
}

int VBfull::main(){
  kmeansRepeat(1);
  
  for(int k=0; k<nMixedElem; k++){
    m_k[k] = mu[k];
  }  
  //    for(int n=0; n<(int)dataTable.size(); n++){
  //      for(int k=0; k<nMixedElem; k++){
  //	//      responsibility[n][k] = 1.0/nMixedElem;
  //	cout << "resp "<<n<<" "<<k<<" "<<responsibility[n][k]<<endl;
  //      }
  //    }

  iteration();
  //cout << "set_ml_parameters"<<endl;
  set_ml_parameters();
  //cout << "get_results_st"<<endl;
  std::cout<<get_results_st();
  
  Write writer(fnOutGaussian);
  writer.open();
  writer.writeLine(get_results_st_tsv());
  writer.close();

  return 0;
}

int VBfull::iteration(){
  double prev_lower_bound = 0;
  int curIteration;

  for(curIteration=0;
      curIteration<=maxIteration; curIteration++){
    //    std::cout << "update_N_xbar_S"<<std::endl;
    update_N_xbar_S();
    //   std::cout << "update_alpha_beta_ny_m_W"<<std::endl;
    update_alpha_beta_ny_m_W();
    //    std::cout << "update_lambda_pi"<<std::endl;
    update_ln_lambda_pi();
    //    std::cout << "update_responsibility"<<std::endl;
    update_responsibility();

    //    std::cout << "cal_lower_bound"<<std::endl;
    double new_lower_bound = cal_lower_bound();
    double diff_lower_bound = new_lower_bound - prev_lower_bound;
    if(printInterval!=0 && curIteration % printInterval == 0){
      std::cout<<"iteration: "<<curIteration<<std::endl;
      std::cout<<"convergence: "<<new_lower_bound<<" "
	       <<diff_lower_bound<<std::endl;
    }
    if(curIteration>0 &&
       fabs(diff_lower_bound)<thresholdConvergence){
      std::cout<<"iteration: "<<curIteration<<std::endl;
      std::cout<<"convergence: "<<new_lower_bound<<" "
	       <<diff_lower_bound<<std::endl;
      break;
    }
    prev_lower_bound=new_lower_bound;
  }
  return 0;

}

int VBfull::set_ml_parameters(){
  for(int k=0;k<nMixedElem;k++){
    pi[k] = (alpha_0 + N_k[k])/(nMixedElem*alpha_0+(double)dataTable.size());
    mu[k] = m_k[k];
    cout<<"ny_k["<<k<<"] "<<ny_k[k]<<" W_k:"<<W_k[k]<<endl;
    invert(ny_k[k] * W_k[k], sigma[k]);
    //    sigma[k] = ny_k[k] * S_k[k];
  }
  return 0;
}

std::string VBfull::get_results_st(){
  std::stringstream ss("");
  for(int k=0;k<nMixedElem;k++){
    cout << "k:"<<k<<endl;
    ss<<"\telement: "<<k<<std::endl;
    ss<<"\t\tpi: "<<pi[k]<<std::endl;
    ss<<"\t\tmu: "<<mu[k]<<std::endl;
    ss<<"\t\tsigma: "<<sigma[k]<<std::endl;
    //    ss<<std::endl;
    //    ss<<"\t\ts:"<<prior_s[k]<<std::endl;
  }
  return ss.str();
}
std::string VBfull::get_results_st_tsv(){
  std::stringstream ss("");
  for(int k=0;k<nMixedElem;k++){
    //cout << k;
    ss<<k<<"\t"<<pi[k];
    for(int l=0; l<dimension; l++){
      ss<<"\t"<<mu[k][l];
    }
    for(int l=0; l<dimension; l++){
      for(int m=0; m<dimension; m++){
	ss<<"\t"<<sigma[k](l,m);
      }
    }
    ss << std::endl;
    //    ss<<std::endl;
    //    ss<<"\t\ts:"<<prior_s[k]<<std::endl;
  }
  return ss.str();
}

