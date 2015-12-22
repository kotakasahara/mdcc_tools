#include "EMAlgorithm.h"
#include "math.h"
namespace ublas=boost::numeric::ublas;

using namespace std;
EMAlgorithm::EMAlgorithm(){}
EMAlgorithm::EMAlgorithm(const ublas::vector<ublas::vector<double> > &inData,
			 std::string inFnOutGaussian,
			 int inNMixedElement,
			 int inMaxIteration,
			 double inThresholdConvergence,
			 int inPrintInterval,
			 int inMaxIterationKmeans){
  dataTable=inData;
  fnOutGaussian = inFnOutGaussian;
  nMixedElem=inNMixedElement;
  maxIteration=inMaxIteration;
  maxIterationKmeans=inMaxIterationKmeans;
  thresholdConvergence=inThresholdConvergence;
  printInterval=inPrintInterval;
  //  maxNKmeans=;
  minDeterminantSigma = 1e-10;
  initialize();
}

int EMAlgorithm::initialize(){
  dimension=dataTable[0].size();
  //  int tmpN=dataTable.size()/nMixedElem;
  //  std::cout << "dimension: "<<dimension<<std::endl;
  //  std::cout << "tmpN: "<<tmpN<<"/"<<dataTable.size()<<std::endl;
  mu=ublas::vector<ublas::vector<double> >(nMixedElem);
  sigma=ublas::vector<ublas::matrix<double> >(nMixedElem);
  pi=ublas::vector<double>(nMixedElem);

  for(int k=0;k<nMixedElem;k++){
    mu[k]=ublas::vector<double>(dimension);
    sigma[k]=ublas::matrix<double>(dimension,dimension);
    pi[k]=0.0;
    for(int j=0;j<dimension;j++){
      mu[k][j]=0.0;
      for(int i=0;i<dimension;i++){
	sigma[k](j,i)=0.0;
      }
    }
  }
  responsibility=ublas::vector<ublas::vector<double> >(dataTable.size());
  for(int n=0;n<(int)dataTable.size();n++){
    responsibility[n]=ublas::vector<double>(nMixedElem);
  }
  
  for(int i=0;i<dimension;i++){
    double min=1e10;
    double max=-1e10;
    //    std::cout<<"init dtsize:"<<dataTable.size()<<std::endl;
    for(int n=0;n<(int)dataTable.size();n++){
      //      std::cout<<"init dt n i:"<<dataTable[n][i]<<std::endl;
      if(dataTable[n][i]>max) max=dataTable[n][i];
      if(dataTable[n][i]<min) min=dataTable[n][i];
    }
    for(int k=0;k<nMixedElem;k++){
      mu[k][i]=((double)rand()/(double)RAND_MAX)*(max-min)+min;
      //      std::cout<<"init mu:"<<mu[k][i]<<" "<<max<<" "<<min<<std::endl;
    }
  }
  
  return 0;
}


int EMAlgorithm::kmeansRepeat(int maxNKmeans){
  std::srand(std::time(NULL));
  bool flgSigmaValid=false;
  int n_kmeans;
  for(n_kmeans=0; n_kmeans<maxNKmeans && !flgSigmaValid; n_kmeans++){
    initialize();
    cout << "kmeans "<<n_kmeans<<endl;
    kmeans();
    flgSigmaValid=true;
    for(int k=0;k<nMixedElem;k++){
      if(determinant(sigma[k])<minDeterminantSigma){
	flgSigmaValid=false;
	break;
      }
    }
  }
  //  if(!flgSigmaValid) exit() ;

  /*
  //calculate initial average
  for(int n=0; n<(int)dataTable.size(); n++){
    //    if(n%tmpN==0){
    //      if(n>0){
    //	for(int j=0; j<dimension;j++){
    //	  mu[k][j]/=(double)(tmpN);
    //	}
    //	k++;
    //	if(k>=nMixedElem) break;
    //      }
    //    }
    for(int k=0;k<nMixedElem;k++){
      for(int j=0;j<dimension;j++){
	mu[k][j]+=dataTable[n][j];
      }
    }
  }
  for(int k=0;k<nMixedElem;k++){
    for(int j=0;j<dimension;j++){
      mu[k][j]/=dataTable.size();
    }
  }

  //calculate initial variance
  for(int n=0; n<(int)dataTable.size(); n++){
    for(int k=0;k<nMixedElem;k++){
      for(int i=0;i<dimension;i++){
	double ei=dataTable[n][i]-mu[k][i];
	for(int j=0;j<dimension;j++){
	  double ej=dataTable[n][j]-mu[k][j];
	  sigma[k](i,j)+=(ei*ej);
	}
      }
    }
  }
  for(int k=0;k<nMixedElem;k++){
    for(int i=0;i<dimension;i++){
      for(int j=0;j<dimension;j++){
	sigma[k](i,j)/=dataTable.size();
      }
    }
  }
  ////////////

  for(int curPi=0;curPi<(int)pi.size();curPi++){
    pi(curPi)=1.0/(double)nMixedElem;
  }
  */

  std::srand(std::time(NULL));
  for(int k=0;k<nMixedElem;k++){
    if(pi[k] == 0 ||
       determinant(sigma[k])<minDeterminantSigma){
      std::cout<<"pi["<<k<<"] set to zero\n";
      for(int i=0;i<dimension;i++){
	for(int j=0;j<dimension;j++){
	  sigma[k](i,j) = 0.0;
	}
	double min=1e10;
	double max=-1e10;
	for(int n=0;n<(int)dataTable.size();n++){
	  if(dataTable[n][i]>max) max=dataTable[n][i];
	  if(dataTable[n][i]<min) min=dataTable[n][i];
	}
	mu[k][i]=((double)rand()/(double)RAND_MAX)*(max-min)+min;
	sigma[k](i,i) = 1.0;
      }
    }
  }
  return n_kmeans;
}

int EMAlgorithm::iteration(){
  //    cout << "pi : "<<pi<<endl;
  double prevLogLikelihood=0;
  int curIteration;
  for(curIteration=0;
      curIteration<=maxIteration; curIteration++){
    ublas::vector<ublas::vector<double> > newResponsibility;
    newResponsibility=stepE();
    responsibility=newResponsibility;
    ublas::vector<ublas::vector<double> > newMu;
    ublas::vector<ublas::matrix<double> > newSigma;
    ublas::vector<double> newPi;
    stepM(newMu, newSigma, newPi);
    mu.assign(newMu);
    sigma.assign(newSigma);
    pi.assign(newPi);
    double newLogLikelihood=calLogLikelihood();
    double diffLogLikelihood=newLogLikelihood-prevLogLikelihood;
    if(curIteration%printInterval==0){
      std::cout<<"iteration: "<<curIteration<<std::endl;
      std::cout<<"convergence: "<<newLogLikelihood<<" "
	       <<diffLogLikelihood<<std::endl;
      std::cout<<outResultsString();
    }
    if(curIteration>0 &&
       fabs(diffLogLikelihood)<thresholdConvergence){
      break;
    }
    prevLogLikelihood=newLogLikelihood;
  }
  std::cout<<"iteration: "<<curIteration<<std::endl;
  std::cout<<outResultsString();

  return curIteration;
}

ublas::vector<ublas::vector<double> > EMAlgorithm::stepE(){
  //   std::cout<<"stepe 1"<<std::endl;
  ublas::vector<ublas::vector<double> > newResponsibility(dataTable.size());
  for(int n=0;n<(int)dataTable.size();n++){
    //    std::cout<<"stepe 2 "<<n<<std::endl;
    newResponsibility[n]=ublas::vector<double>(nMixedElem);
    ublas::vector<double> num(nMixedElem);
    double den=0;
    for(int j=0;j<nMixedElem;j++){
      //      std::cout << "stepe 2.1 j:"<<pi[j]<<std::endl;
      //      std::cout << "stepe 2.1.2 sigme "<<sigma[j]<<endl;
      num[j]=pi[j]*calGaussian(dataTable[n],mu[j],sigma[j]);
      //      std::cout << "stepe 2.2 j:"<<num[j]<<std::endl;
      den+=num[j];
    }
    //    std::cout<<"stepe 3 n="<<n<<" den="<<den<<std::endl;
    for(int k=0;k<nMixedElem;k++){
      newResponsibility[n][k]=num[k]/den;
      //      std::cout<<"stepe 4 gamma["<<n<<"]["<<k<<"]="<<newGamma[n][k]<<std::endl;
    }
  }
  return newResponsibility;
}

int EMAlgorithm::stepM(ublas::vector<ublas::vector<double> > &newMu,
		       ublas::vector<ublas::matrix<double> > &newSigma,
		       ublas::vector<double> &newPi){  
  std::srand(std::time(NULL));

  newMu=ublas::vector<ublas::vector<double> >(nMixedElem);
  newSigma=ublas::vector<ublas::matrix<double> >(nMixedElem);
  newPi=ublas::vector<double>(nMixedElem);
  double sumNk=0;
  for(int k=0;k<nMixedElem;k++){
    newMu[k]=ublas::vector<double>(dimension);
    newMu[k]*=0;
    double Nk=0;
    for(int n=0;n<(int)dataTable.size();n++){
      Nk+=responsibility[n][k];
      newMu[k]+=responsibility[n][k]*dataTable[n];
    }
    sumNk+=Nk;
    //    cout << "m0 mu["<<k<<"]: "<<newMu<<endl;
    newMu[k]/=Nk;
    newSigma[k]=ublas::matrix<double>(dimension,dimension);
    newSigma[k]*=0;
    //    std::cout << "m1 newSigma "<<newSigma.size()<<" "<<newSigma[0].size1()<<" "<<newSigma[0].size2()<<std::endl;
    //    cout <<"stepM newSigma[k]:" <<newSigma[k]<<endl;
    for(int n=0;n<(int)dataTable.size();n++){
      ublas::vector<double> hensaV=dataTable[n]-newMu[k];
      //      cout<<"stepM n:"<<n<<" hensaV:"<<hensaV<<endl;
      ublas::matrix<double> hensa=vectorToMatrix(hensaV);
      //      cout<<"stepM n:"<<n<<" hensa:"<<hensa<<endl;
      ublas::matrix<double> hensaT=vectorToMatrixTranspose(hensaV);
      //      cout<<"stepM n:"<<n<<" hensaT:"<<hensaT<<endl;
      newSigma[k]+=responsibility[n][k]*prod(hensaT,hensa);
      //      cout <<"stepM newSigma[k]:" <<newSigma[k]<<endl;
      //      cout <<"stepM responsibility[n][k]:"<<responsibility[n][k]<<endl;
      //      cout <<"stepM prod :"<<prod(hensaT,hensa)<<endl;
      //      cout <<"stepM responsibility[n][k]*prod:"<<responsibility[n][k]*prod(hensaT,hensa)<<endl;
    }
    newSigma[k]/=Nk;
    //    std::cout << "m2 newSigma "<<newSigma.size()<<" "<<newSigma[0].size1()<<" "<<newSigma[0].size2()<<std::endl;
    //    std::cout << newSigma <<endl;
    newPi[k]=(double)Nk/(double)dataTable.size();
    while(determinant(newSigma[k])<minDeterminantSigma){
      for(int i=0;i<dimension;i++){
	double min=1e10;
	double max=-1e10;
	for(int n=0;n<(int)dataTable.size();n++){
	  if(dataTable[n][i]>max) max=dataTable[n][i];
	  if(dataTable[n][i]<min) min=dataTable[n][i];
	}
	newMu[k][i]=((double)rand()/(double)RAND_MAX)*(max-min)+min;
      }
      newSigma[k] *= 0;
      for(int n=0;n<(int)dataTable.size();n++){
	ublas::vector<double> hensaV=dataTable[n]-newMu[k];
	ublas::matrix<double> hensa=vectorToMatrix(hensaV);
	ublas::matrix<double> hensaT=vectorToMatrixTranspose(hensaV);
	newSigma[k]+=prod(hensaT,hensa);
      }
    }
  }
  return 0;
}

double EMAlgorithm::calLogLikelihood(){
  double logLikelihood=0;
  //  std::cout<<"likelihood"<<std::endl;
  for(int n=0;n<(int)dataTable.size();n++){
    double tmpLikelihood=0;
    for(int k=0;k<nMixedElem;k++){
      //      std::cout<<"likelihood n="<<n<<" k="<<k<<std::endl;
      tmpLikelihood+=pi[k]*calGaussian(dataTable[n],mu[k],sigma[k]);
      //      cout << "like : pi : "<<pi[k]<<endl;
      //      cout << "like : mu : "<<mu[k]<<endl;
      //      cout << "like : sigma : "<<sigma[k]<<endl;
      //      cout << "like : g  : "<<calGaussian(dataTable[n],mu[k],sigma[k])<<endl;
    }
    //    std::cout<<"loglike "<<tmpLikelihood<<std::endl;
    logLikelihood+=log(tmpLikelihood);
  }
  return logLikelihood;
}


void EMAlgorithm::test(){
  //  lu_substitute(gSigma,pm,gSigmaInv);
  /*  
 typedef boost::numeric::ublas::ublas::vector<double> dublas::vector;
 typedef boost::numeric::ublas::ublas::matrix<double> dublas::matrix;
 dublas::matrix B = identity_ublas::matrix<double>(3);
 dublas::matrix A(3,3);
 A(0,0) = 1.0;  A(0,1) = 2.0;  A(0,2) = 3.0;
 A(1,0) = 2.0;  A(1,1) = 1.0;  A(1,2) = -3.0;
 A(2,0) = 4.0;  A(2,1) = -3.0; A(2,2) = 2.0;
 dublas::matrix CopyOfA(A);
 permutation_ublas::matrix<dublas::matrix::size_type> pm(A.size1());
 std::cout << "A = " << A << std::endl;
 std::cout << "B = " << B << std::endl;
 lu_factorize(A,pm);
 lu_substitute(A,pm,B);
 std::cout.precision(4);
 std::cout << "Ainv X = " << B << std::endl;
 std::cout << "A * X = " << prod(CopyOfA,B) << std::endl;
  
  std::cout<<BOOST_LIB_VERSION<<std::endl;
  */
}

std::string EMAlgorithm::outResultsString(){
  std::stringstream ss;
  for(int k=0;k<nMixedElem;k++){
    ss<<"\telement: "<<k<<std::endl;
    ss<<"\t\tpi: "<<pi[k]<<std::endl;
    ss<<"\t\tmu: "<<mu[k]<<std::endl;
    ss<<"\t\tsigma: "<<sigma[k]<<std::endl;
    //    ss<<std::endl;
  }
  double mdl=calMDL();
  ss<<"mdl: "<<mdl<<endl;
  return ss.str();
}

double EMAlgorithm::calMDL(){
  double mdl_1=0;
  for(int n=0;n<(int)dataTable.size();n++){
    double prob=0;
    for(int k=0;k<nMixedElem;k++){
      prob += pi[k]*calGaussian(dataTable[n],mu[k],sigma[k]);
      //      cout << "mdl_1:n="<<n<<" k="<<k<<" : ";
      //      cout << prob<<endl;
    }
    mdl_1 += log(prob);
  }
  mdl_1 *= -2;
  double mdl_2 = nMixedElem*(dimension+1)*(dimension+2) * 0.5 - 1;
  //  cout << "mdl_1:"<<mdl_1<<" mdl_2:"<<mdl_2<<endl;
  mdl_2 *= log(dataTable.size());
  return mdl_1+mdl_2;
}

int EMAlgorithm::kmeans(){
   std::cout << "kmeans"<<std::endl;
  int curIteration;
  bool flgChange=false;
  ublas::vector<int> nMember(nMixedElem);
  ublas::vector<ublas::vector<double> > prevResponsibility(responsibility);
  for(curIteration=0;curIteration<maxIterationKmeans;curIteration++){
    //        std::cout << "kmeans it:"<<curIteration<<std::endl;
    for(int n=0;n<(int)dataTable.size();n++){
      double minD=1e10;
      for(int k=0;k<nMixedElem;k++){
	ublas::vector<double> hensaV=dataTable[n]-mu[k];
	ublas::matrix<double> hensa=vectorToMatrix(hensaV);
	ublas::matrix<double> hensaT=vectorToMatrixTranspose(hensaV);
	double curD;
	ublas::matrix<double> tmp=prod(hensa,hensaT);
	curD=tmp(0,0);
	if(curD<minD){
	  minD=curD;
	  for(int j=0;j<nMixedElem;j++){
	    if(k==j) responsibility[n][j]=1.0;
	    else responsibility[n][j]=0.0;
	  }
	}
      }
    }
    for(int k=0;k<nMixedElem;k++){
      mu[k]=ublas::vector<double>(dimension);
      nMember[k]=0;
    }
    for(int n=0;n<(int)dataTable.size();n++){
      for(int k=0;k<nMixedElem;k++){
	if(responsibility[n][k]>0){
	  nMember[k]++;
	  mu[k]+=dataTable[n];
	}
      }
    }
    for(int k=0;k<nMixedElem;k++){
      mu[k]/=nMember[k];
    }
    flgChange=false;
    //    for(int n=0;n<(int)dataTable.size()&&!flgChange;n++){
    //      for(int k=0;k<nMixedElem&&!flgChange;k++){
    //	if(responsibility[n][k]!=prevResponsibility[n][k]) flgChange=true;
    //      }
    //    }
    //    if(!flgChange){break;}
    //    prevResponsibility=responsibility;
  }

  for(int k=0;k<nMixedElem;k++){
    for(int n=0;n<(int)dataTable.size();n++){
      if(responsibility[n][k]==1){
	for(int i=0;i<dimension;i++){
	  double ei=dataTable[n][i]-mu[k][i];
	  for(int j=0;j<dimension;j++){
	    double ej=dataTable[n][j]-mu[k][j];
	    sigma[k](i,j)+=ei*ej;
	  }
	}
      }
    }
  }
  for(int k=0;k<nMixedElem;k++){
    sigma[k]/=nMember[k];
    pi[k]=(double)nMember[k]/(double)dataTable.size();
  }
  cout <<"nMember:"<< nMember<<endl;;
  cout <<"pi:"<< pi<<endl;
  cout <<"mu:"<< mu<<endl;
  cout <<"sigma:"<<sigma<<endl;
  return 0;

}

