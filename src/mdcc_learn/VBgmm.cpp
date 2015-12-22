#include "VBgmm.h"
#include "math.h"
VBgmm::VBgmm() : EMAlgorithm(){
}

VBgmm::VBgmm(const ublas::vector<ublas::vector<double> > &inData,
	     std::string inFnOutGaussian,
	     int inNMixedElement,
	     int inMaxIteration,
	     double inThresholdConvergence,
	     int inPrintInterval,
	     double in_hp_beta,
	     int inMaxIterationKmeans)
  : EMAlgorithm(inData, inFnOutGaussian,
		inNMixedElement, inMaxIteration,
		inThresholdConvergence, inPrintInterval,
		inMaxIterationKmeans){
  hp_beta = in_hp_beta;
}

int VBgmm::initialize(){
  prior_s = ublas::vector<ublas::vector<double> > (nMixedElem);
  mean_prior_mean = ublas::vector<ublas::vector<double> >(nMixedElem);
  prec_prior_mean = ublas::vector<ublas::matrix<double> >(nMixedElem);
  dof_prior_prec = ublas::vector<double>(nMixedElem);
  scale_prior_prec = ublas::vector<ublas::matrix<double> >(nMixedElem);

  ublas::matrix<double> unitmatrix(dimension,dimension);
  for(int d1=0;d1<dimension;d1++){
    for(int d2=0;d2<dimension;d2++){
      unitmatrix(d1,d2) = 0.0;
    }
    unitmatrix(d1,d1) = 1.0;
  }
  for(int i=0;i<nMixedElem;i++){
    mean_prior_mean[i] = mu[i];
    prec_prior_mean[i] = hp_beta * unitmatrix;
    dof_prior_prec[i] = (int)(dimension*(dimension+1.0)/2.0);
    scale_prior_prec[i] = hp_beta * unitmatrix;
    pi[i] = 1.0/nMixedElem;
    //    std::cout <<"initial pi["<<i<<"] "<<pi[i]<<std::endl;
    //    std::cout <<"init mean_prior_mean["<<i<<"] "
    //	      <<mean_prior_mean[i]<<std::endl;
  }
  //initialize  prior_s
  double tmp_prior_s = 1.0/nMixedElem;
  for(int i=0;i<nMixedElem;i++){
    prior_s[i] = ublas::vector<double>(dataTable.size());
    for(int n=0;n<(int)dataTable.size();n++){
      prior_s[i][n] = tmp_prior_s;
    }
  }
  hp_prec_g = hp_beta * unitmatrix;
  //  hp_dof_w = (int)(dimension*(dimension+1.0)/2.0);
  hp_dof_w = (int)(dimension+2);
  hp_scale_w = hp_beta * unitmatrix;

  return 0;
}

int VBgmm::main(){
  std::cout << "kmeans"<<std::endl;
  kmeansRepeat(1);
  initialize();
  //  for(int i=0;i<nMixedElem;i++){
  //  std::cout <<"mu["<<i<<"] "
  //	    <<mu[i]<<std::endl;
  //  }
  iteration();
  set_ml_parameters();
  std::cout<<get_results_st();
  return 0;
}

int VBgmm::iteration(){
  double prev_lower_bound = 0;
  int curIteration;

  for(curIteration=0;
      curIteration<=maxIteration; curIteration++){
    //    std::cout << "step e"<<std::endl;
    stepE();
    //    std::cout << "step m"<<std::endl;
    stepM();
    //    std::cout << "cal_lower_bound"<<std::endl;
    double new_lower_bound = cal_lower_bound();
    double diff_lower_bound = new_lower_bound - prev_lower_bound;
    if(curIteration % printInterval == 0){
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
int VBgmm::stepE(){
  stepE_prior_for_prec();
  stepE_prior_for_mean();
  return 0;
}

int VBgmm::stepE_prior_for_mean(){
  for(int i=0;i<nMixedElem;i++){  
    if(pi[i] <= 1e-5) continue;
    //ex_prec <Ti>
    //  std::cout<<"stepE_prior_for_mean  i="<<i<<std::endl;
    //  std::cout<<"ex_prec"<<std::endl;
    ublas::matrix<double> ex_prec(dimension,dimension);
    ublas::matrix<double> inv_scale(dimension,dimension);
    invert(scale_prior_prec[i], inv_scale);
    ex_prec = dof_prior_prec[i]*inv_scale;
    //sum_prior_s;
    //    std::cout<<"sum_prior_s"<<std::endl;
    double sum_prior_s = 0;
    ublas::vector<double> sum_x_prior_s(dimension);
    for(int s=0;s<dimension;s++) sum_x_prior_s[s]=0.0;
    for(int n=0;n<(int)dataTable.size();n++){
      //       std::cout<<dataTable[n].size()<<std::endl;
      sum_prior_s += prior_s[i][n];
      //      std::cout<<sum_prior_s<<std::endl;
      sum_x_prior_s += dataTable[n] * prior_s[i][n];
      //      std::cout<<sum_x_prior_s<<std::endl;
    }
    //    std::cout<<hp_prec_g.size1()<<" "
    //	     <<hp_prec_g.size2()<<" "
    //	     <<ex_prec.size1()<<" "
    //	     <<ex_prec.size2()<<std::endl;
    prec_prior_mean[i] = hp_prec_g + ex_prec * sum_prior_s;
    
   //;
    //    std::cout<<"sum_x_prior_s"<<std::endl;
    ublas::matrix<double> prec_inv(dimension,dimension);
    invert(prec_prior_mean[i], prec_inv);
    ublas::vector<double> ex_prec_sum_x_prior_s
      = ublas::prod(ex_prec,sum_x_prior_s);
    
    //    std::cout<<"mean_prior_mean"<<std::endl;
    mean_prior_mean[i] = ublas::prod(prec_inv, ex_prec_sum_x_prior_s);
    //    std::cout<<"stepE_prior_for_mean  i="<<i<<" "
    //	     <<mean_prior_mean[i]<<std::endl;
  }
  return 0;
}

int VBgmm::stepE_prior_for_prec(){
  for(int i=0;i<nMixedElem;i++){  
    if(pi[i] <= 1e-5) continue;
    //    std::cout<<"stepE_prior_for_prec  i="<<i<<std::endl;
    // define and initialize
    double sum_prior_s = 0;
    ublas::matrix<double> sum_xxs(dimension,dimension);
    ublas::matrix<double> sum_xsmu(dimension,dimension);
    ublas::matrix<double> sum_xs(1,dimension);
    for(int a=0;a<dimension;a++){
      sum_xs(0,a) = .0;
      for(int b=0;b<dimension;b++){
	sum_xxs(a,b) = .0;
	sum_xsmu(a,b) = .0;
      }
    }
    // cal summation among data points
    //    std::cout<<"sum"<<std::endl;
    ublas::matrix<double> ex_mean
      = vectorToMatrix(mean_prior_mean[i]);
    ublas::matrix<double> ex_mean_t
      = vectorToMatrixTranspose(mean_prior_mean[i]);
    ublas::matrix<double> prec_inv(dimension,dimension);
    invert(prec_prior_mean[i],prec_inv);
    //    ublas::matrix<double> prd = ublas::prod(ex_mean,ex_mean_t);
    //    std::cout<<"dbg d "<<prd.size1()<<" "<<prd.size2()<<std::endl;
    ublas::matrix<double> ex_mean_2
      = prec_inv + ublas::prod(ex_mean_t,ex_mean);
    for(int n=0;n<(int)dataTable.size();n++){
      sum_prior_s += prior_s[i][n];
      ublas::matrix<double> data_n = vectorToMatrix(dataTable[n]);
      ublas::matrix<double> data_n_t = vectorToMatrixTranspose(dataTable[n]);
      sum_xxs += ublas::prod(data_n_t,data_n) * prior_s[i][n];
      sum_xsmu += ublas::prod(data_n_t*prior_s[i][n], ex_mean);
      sum_xs += data_n * prior_s[i][n];
    }
    
    //    std::cout<<"dof, scale"<<std::endl;
    dof_prior_prec[i] = hp_dof_w + sum_prior_s;
    scale_prior_prec[i] = hp_scale_w
      + sum_xxs
      - sum_xsmu
      - ublas::prod(ex_mean_t,sum_xs)
      + ex_mean_2*sum_prior_s;
  }
  return 0;
}

int VBgmm::stepM(){
  //update prior_s (p) 
  ublas::vector<ublas::vector<double> > prior_s_sub(nMixedElem);
  ublas::vector<double> sum_prior_s_sub(dataTable.size());
  
  for(int i=0;i<nMixedElem;i++){  
    if(pi[i] <= 1e-5) continue;
    prior_s_sub[i] = ublas::vector<double>(dataTable.size());
    //ex_ln_prec <ln|Ti|>
    double ex_ln_prec_1 = 0.0;
    for(int s=1;s<=dimension;s++){
      ex_ln_prec_1 += digamma((dof_prior_prec[i]+1.0-s)/2.0);
    }
    double ex_ln_prec = ex_ln_prec_1 + dimension*log(2.0)
      - log(determinant(scale_prior_prec[i]));
    //ex_prec <Ti>
    ublas::matrix<double> ex_prec(dimension,dimension);
    ublas::matrix<double> inv_scale(dimension,dimension);
    invert(scale_prior_prec[i], inv_scale);
    ex_prec = dof_prior_prec[i]*inv_scale;
    //<mu i>
    //<prod<mu,muT>>
    ublas::matrix<double> prec_inv(dimension,dimension);
    invert(prec_prior_mean[i],prec_inv);
    ublas::matrix<double> ex_mean
      = vectorToMatrix(mean_prior_mean[i]);
    ublas::matrix<double> ex_mean_t
      = vectorToMatrixTranspose(mean_prior_mean[i]);
    ublas::matrix<double> ex_mean_2
      = prec_inv + ublas::prod(ex_mean_t,ex_mean);
    
    for(int n=0;n<(int)dataTable.size();n++){
      ublas::matrix<double> data_n
	= vectorToMatrix(dataTable[n]);
      ublas::matrix<double> data_n_t
	= vectorToMatrixTranspose(dataTable[n]);
      //prod<xn,xnT>
      ublas::matrix<double> xnxn = ublas::prod(data_n_t,data_n);
      //prod<<mu>,xnT>
      ublas::matrix<double> muxn = ublas::prod(ex_mean_t,data_n);
      //prod<xn,muT>
      ublas::matrix<double> xnmu = ublas::prod(data_n_t,ex_mean);

      //      std::cout<<"uhoho "
      //	       <<ex_ln_prec<<" "
      //	       <<pi[i]<<" "
      //	       <<log(pi[i])<<" "
      //	       <<trace(ublas::prod(ex_prec,(xnxn-muxn-xnmu+ex_mean_2)))
      //	       <<std::endl;
      prior_s_sub[i][n]
	= exp( 0.5 * ex_ln_prec
	       + log(pi[i])
	      - 0.5*trace(ublas::prod(ex_prec,(xnxn-muxn-xnmu+ex_mean_2))));
      //      sum_prior_s_sub[i] += prior_s_sub[i][n];
    }
  }
  for(int n=0;n<(int)dataTable.size();n++){
    sum_prior_s_sub[n] = 0.0;
    for(int i=0;i<nMixedElem;i++){  
    if(pi[i] <= 1e-5) continue;
      sum_prior_s_sub[n] += prior_s_sub[i][n];
    }
  }

  
  //  double tmp=0;
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    for(int n=0;n<(int)dataTable.size(); n++){
      //      std::cout<<"pri_s update "
      //	       <<" i:"<<i<<" n:"<<n<<" "
      //	       <<prior_s_sub[i][n]<<" "
      //	       <<sum_prior_s_sub(i)<<" "
      //	       <<std::endl;
      prior_s[i][n] = prior_s_sub[i][n] / sum_prior_s_sub(n);
      //      tmp += prior_s[i][n];
    }
  }
  //  std::cout<<"sum prior_s = "<<tmp<<std::endl;
  
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    pi[i] = 0.0;
    for(int n=0;n<(int)dataTable.size(); n++){
      pi[i] += prior_s[i][n];
    }
    pi[i] /= (double)dataTable.size();
  }
  return 0;
}

double VBgmm::cal_lower_bound(){

  ublas::vector<ublas::matrix<double> > ex_prec(nMixedElem);
  ublas::vector<ublas::matrix<double> > inv_scale(nMixedElem);
  ublas::vector<double> ex_ln_prec(nMixedElem);
  ublas::vector<ublas::matrix<double> > ex_mean_2(nMixedElem);

  int curNMixedElem=0;
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] > 1e-5) curNMixedElem++;
  }
  
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    //ex_prec <Ti>
    invert(scale_prior_prec[i], inv_scale[i]);
    ex_prec[i] = dof_prior_prec[i]*inv_scale[i];

    //ex_ln_prec <ln|Ti|>
    double ex_ln_prec_1 = 0;
    for(int d=1;d<=dimension;d++){
      ex_ln_prec_1 += digamma(dof_prior_prec[i]+1.0-d)/2.0;
    }
    ex_ln_prec[i] = ex_ln_prec_1 + dimension*log(2.0)
      - log(determinant(scale_prior_prec[i]));
    
    //ex_mumu
    ublas::matrix<double> prec_inv(dimension,dimension);
    invert(prec_prior_mean[i], prec_inv);
    ublas::matrix<double> ex_mean = vectorToMatrix(mean_prior_mean[i]);
    ublas::matrix<double> ex_mean_t
      = vectorToMatrixTranspose(mean_prior_mean[i]);
    ex_mean_2[i] = prec_inv + ublas::prod(ex_mean_t,ex_mean);
  }
  
  //eq 28
  double ex_ln_p_d_given_mu_t_s = 0;
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    ublas::matrix<double> ex_mean = vectorToMatrix(mean_prior_mean[i]);
    ublas::matrix<double> ex_mean_t
      = vectorToMatrixTranspose(mean_prior_mean[i]);
    for(int n=0;n<(int)dataTable.size();n++){
      ublas::matrix<double> data_n
	= vectorToMatrix(dataTable[n]);
      ublas::matrix<double> data_n_t
	= vectorToMatrixTranspose(dataTable[n]);
      ublas::matrix<double> tmp
	= ublas::prod(data_n_t,data_n)
	- ublas::prod(data_n_t,ex_mean)
	- ublas::prod(ex_mean_t, data_n)
	+ ex_mean_2[i];
      ublas::matrix<double> tmp2 = ublas::prod(ex_prec[i],tmp);
      //      std::cout<<"pri_s "<<prior_s[i][n]
      //	       <<" "<<ex_ln_prec[i]
      //	       <<" "<<trace(tmp2)<<std::endl;
      ex_ln_p_d_given_mu_t_s += prior_s[i][n]
	*(0.5*ex_ln_prec[i]
	  - dimension*0.5*log(2*PI)
	  - 0.5*trace(tmp2));
    }
  }
  //eq 29
  double ex_ln_p_s = 0;
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    for(int n=0;n<(int)dataTable.size();n++){
      ex_ln_p_s += prior_s[i][n]*log(pi[i]);
    }
  }

  //eq 30
  double ex_ln_p_mu = curNMixedElem*dimension*log(hp_beta/(2*PI));
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    ublas::matrix<double> ex_mean
      = vectorToMatrix(mean_prior_mean[i]);
    ublas::matrix<double> ex_mean_t
      = vectorToMatrixTranspose(mean_prior_mean[i]);
    double tmp_mumu = determinant(inv_scale[i])+
      ublas::prod(ex_mean,ex_mean_t)(0,0);

    ex_ln_p_mu -= (hp_beta*0.5)*tmp_mumu;
  }
  
  //eq 31
  double ex_ln_p_t=0;
  double sum_ex_ln_t=0;
  ublas::matrix<double> sum_ex_t(dimension,dimension);
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    sum_ex_ln_t += ex_ln_prec[i];
    sum_ex_t += ex_prec[i];
  }
  double tmp_ln_gam=0;
  for(int s=0;s<dimension;s++){
       tmp_ln_gam += lgamma((hp_dof_w+1.0-s)*0.5);
    //    tmp_ln_gam += log( gamma(hp_dof_w+1.0-s)*0.5 );
    //    tmp_ln_gam += log(  );
  }

  ex_ln_p_t = curNMixedElem
    * (-0.5*hp_dof_w*dimension*log(2.0)
       -dimension*(dimension-1.0)*0.25*log(PI)
       -tmp_ln_gam
       +0.5*hp_dof_w*log(determinant(hp_scale_w)))
    + 0.5*(hp_dof_w-dimension-1.0)*sum_ex_ln_t
    - 0.5*trace(ublas::prod(hp_scale_w,sum_ex_t));
    
  //eq 32
  double ex_ln_q_s = 0;
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    for(int n=0;n<(int)dataTable.size();n++){
		std::cout<<"["<<i<<"]["<<n<<"] pri_s="<<prior_s[i][n]<<std::endl;
		ex_ln_q_s += prior_s[i][n]*log(prior_s[i][n]);
    }
  }

  //eq 33
  double ex_ln_q_mu = 0;
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    ex_ln_q_mu += -0.5*dimension*(1.0+log(2.0*PI))
      + 0.5*log(determinant(prec_prior_mean[i]));
  }
  
  //eq 34
  double ex_ln_q_t = 0;
  double term2 = -0.25*dimension*(dimension-1.0)*log(PI);
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    double term1 = -0.5*dof_prior_prec[i]*dimension*log(2.0);
    double term3 = 0;
    for(int s=0;s<dimension;s++){
      term3 += lgamma(0.5*(dof_prior_prec[i]+1.0-s));
    }
    double term4 = 0.5*dof_prior_prec[i]*log(determinant(scale_prior_prec[i]));
    double term5 = 0.5*(dof_prior_prec[i]-dimension-1.0)*ex_ln_prec[i];
    double term6 = -0.5*trace(ublas::prod(scale_prior_prec[i],ex_prec[i]));
    ex_ln_q_t += term1 + term2 + term3 + term4 + term5 + term6;
  }
    std::cout<<"cal_lower "<<ex_ln_p_d_given_mu_t_s
  	   <<" " << ex_ln_p_s
     <<" " <<     ex_ln_p_mu
     <<" " <<     ex_ln_p_t
     <<" " <<     ex_ln_q_s
     <<" " <<     ex_ln_q_mu
     <<" " <<     ex_ln_q_t
  	   <<std::endl;
  
  return ex_ln_p_d_given_mu_t_s
    + ex_ln_p_s
    + ex_ln_p_mu
    + ex_ln_p_t
    - ex_ln_q_s
    - ex_ln_q_mu
    - ex_ln_q_t;
}

int VBgmm::set_ml_parameters(){
  for(int i=0;i<nMixedElem;i++){
    if(pi[i] <= 1e-5) continue;
    mu[i] = mean_prior_mean[i];
    invert(dof_prior_prec[i] * scale_prior_prec[i], sigma[i]);
  }
  return 0;
}

std::string VBgmm::get_results_st(){
  std::stringstream ss;
  for(int k=0;k<nMixedElem;k++){
    ss<<"\telement: "<<k<<std::endl;
    ss<<"\t\tpi: "<<pi[k]<<std::endl;
    ss<<"\t\tmu: "<<mu[k]<<std::endl;
    ss<<"\t\tsigma: "<<sigma[k]<<std::endl;
    //    ss<<std::endl;
    //    ss<<"\t\ts:"<<prior_s[k]<<std::endl;
  }
  for(int k=0;k<nMixedElem;k++){
    ss<<"element: "<<k<<std::endl;
    ss<<"mean_prior_mean: "<<mean_prior_mean[k]<<std::endl;
    ss<<"prec_prior_mean:"<<prec_prior_mean[k]<<std::endl;
    ss<<"dof_prior_prec:"<<dof_prior_prec[k]<<std::endl;
    ss<<"scale_prior_prec:"<<scale_prior_prec[k]<<std::endl;
  }
  return ss.str();
}
