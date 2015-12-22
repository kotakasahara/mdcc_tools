#include "Gaussian.h"
#include "math.h"


/////////////////// Gaussian /////////////////////

Gaussian::Gaussian(){
}
Gaussian::Gaussian(int inId, string inType,
		   ublas::vector<double> inMu, ublas::matrix<double> inSigma){
  id = inId;
  type = inType;
  mu = inMu;
  sigma = inSigma;
  invert(sigma, sigma_inv);
  sigma_det = determinant(sigma);
}
double Gaussian::calProb(ublas::vector<double> x){
  return calGaussianInv(x, mu, sigma_inv, sigma_det);
}

double Gaussian::setProbabilityMesh(ublas::vector<double> init,
				    ublas::vector<double> stepwidth,
				    ublas::vector<int> n_steps,
				    bool xyzToSpherical,
				    double laplace_delta){
  //  cout << "GaussianMixture::setProbabilityMesh"<<endl;
  int dimension = init.size();
  ublas::vector<double> var = init;
  ublas::vector<int> cur_step(dimension);
  for(int i=0; i<dimension; i++) cur_step[i]=0;
  int sum_n_steps=1;
  ublas::vector<int>::iterator itr_nsteps;
  for(itr_nsteps=n_steps.begin();
      itr_nsteps!=n_steps.end(); itr_nsteps++){
    sum_n_steps *= *itr_nsteps;
  }
  if (sum_n_steps == 0) {cout <<"sum_n_steps = 0!!"<<endl; exit(1);}
  //  cout << "sum_n_steps "<<sum_n_steps<<endl;
  mesh = ublas::vector<double>(sum_n_steps);
  for(int i=0;i<sum_n_steps;i++){mesh[i]=0.0;}
  mesh_sum = 0.0;
  for(int i=0; i<sum_n_steps; i++){
    int mlt = 1;
    for(int d=0; d<n_steps.size();d++){
      cur_step[d] = (i/mlt)%n_steps[d];
      var[d] = init[d]+cur_step[d]*stepwidth[d];
      mlt *= n_steps[d];
    }
    if(xyzToSpherical){
      Coord var_c(var[0],var[1],var[2]);
      if(var[0]==0 || var[1]==0 || var[2]==0) {
	continue;
      }
      var[0] = var_c.distTo(Coord(0,0,0));
      var[1] = var_c.calIgTheta();
      var[2] = var_c.calIgPhi();
    }
    double prob = calProb(var);
    //    if(prob < 1e-10){
      //      cout << "SMALL " << prob <<endl;
    //    }
    prob += laplace_delta;
    mesh[i] = prob;
    mesh_sum += prob;
  }
  //  cout << "setProbabilityMesh sum_mesh "<<sum_mesh <<endl;
  for(int i=0;i<sum_n_steps;i++){
    mesh[i] /= mesh_sum;
  }
  return 0;
}
double Gaussian::setProbabilityMeshSum(ublas::vector<double> init,
				   ublas::vector<double> stepwidth,
				   ublas::vector<int> n_steps,
				   bool xyzToSpherical,
				   double laplace_delta){
  //  cout << "GaussianMixture::setProbabilityMesh"<<endl;
  int dimension = init.size();
  ublas::vector<double> var = init;
  ublas::vector<int> cur_step(dimension);
  for(int i=0; i<dimension; i++) cur_step[i]=0;
  int sum_n_steps=1;
  ublas::vector<int>::iterator itr_nsteps;
  for(itr_nsteps=n_steps.begin();
      itr_nsteps!=n_steps.end(); itr_nsteps++){
    sum_n_steps *= *itr_nsteps;
  }
  if (sum_n_steps == 0) {cout <<"sum_n_steps = 0!!"<<endl; exit(1);}
  //  cout << "sum_n_steps "<<sum_n_steps<<endl;
  mesh_sum = 0.0;
  for(int i=0; i<sum_n_steps; i++){
    int mlt = 1;
    for(int d=0; d<n_steps.size();d++){
      cur_step[d] = (i/mlt)%n_steps[d];
      var[d] = init[d]+cur_step[d]*stepwidth[d];
      mlt *= n_steps[d];
    }
    if(xyzToSpherical){
      Coord var_c(var[0],var[1],var[2]);
      if(var[0]==0 || var[1]==0 || var[2]==0) {
	continue;
      }
      var[0] = var_c.distTo(Coord(0,0,0));
      var[1] = var_c.calIgTheta();
      var[2] = var_c.calIgPhi();
    }
    double prob = calProb(var);
    //    if(prob < 1e-10){
      //      cout << "SMALL " << prob <<endl;
    //    }
    prob += laplace_delta;
    mesh_sum += prob;
  }
  return mesh_sum;
}
double Gaussian::calGaussMaharanobisDist(ublas::vector<double> x){
  return calMaharanobisDistInv(x, mu, sigma_inv);
}

double Gaussian::calKLDivergenceMesh(Gaussian *ge){
  //  cout << "GaussianMixture::calKLDivergenceMesh"<<endl;

  double kld = 0.0;
  for(int i=0;i<(int)mesh.size();i++){
    double p1n = mesh[i];
    double p2n = ge->getMesh(i);
    if(p1n>1e-15 && p2n>1e-15)
      kld += p1n*log(p1n/p2n);
    else
      cerr<< "probability is too small. "<<i<<" "<<p1n<<" "<<p2n
	  << " " << id << " " <<ge->getId()<<endl;
  }
  return kld;
}
double Gaussian::calKlDivergence(ublas::vector<double> init,
				 ublas::vector<double> stepwidth,
				 ublas::vector<int> n_steps,
				 bool xyzToSpherical,
				 double laplace_delta,
				 Gaussian* ge){
  //  cout << "GaussianMixture::setProbabilityMesh"<<endl;
  int dimension = init.size();
  ublas::vector<double> var = init;
  ublas::vector<int> cur_step(dimension);
  for(int i=0; i<dimension; i++) cur_step[i]=0;
  int sum_n_steps=1;
  ublas::vector<int>::iterator itr_nsteps;
  for(itr_nsteps=n_steps.begin();
      itr_nsteps!=n_steps.end(); itr_nsteps++){
    sum_n_steps *= *itr_nsteps;
  }
  if (sum_n_steps == 0) {cout <<"sum_n_steps = 0!!"<<endl; exit(1);}
  //  cout << "sum_n_steps "<<sum_n_steps<<endl;
  double kld = 0.0;
  for(int i=0; i<sum_n_steps; i++){
    int mlt = 1;
    for(int d=0; d<n_steps.size();d++){
      cur_step[d] = (i/mlt)%n_steps[d];
      var[d] = init[d]+cur_step[d]*stepwidth[d];
      mlt *= n_steps[d];
    }
    if(xyzToSpherical){
      Coord var_c(var[0],var[1],var[2]);
      if(var[0]==0 || var[1]==0 || var[2]==0) {
	continue;
      }
      var[0] = var_c.distTo(Coord(0,0,0));
      var[1] = var_c.calIgTheta();
      var[2] = var_c.calIgPhi();
    }
    double prob1 = (calProb(var)+laplace_delta) / mesh_sum;
    double prob2 = (ge->calProb(var)+laplace_delta) / ge->getMeshSum();
    if(prob1>1e-15 && prob2>1e-15)
      kld += prob1*log(prob1/prob2);
    else
      cerr<< "probability is too small. "<<i<<" "<<prob1<<" "<<prob2
	  << " " << id << " " <<ge->getId()<<endl;
    //    if(prob < 1e-10){
      //      cout << "SMALL " << prob <<endl;
    //    }
  }
  return kld;
}
double Gaussian::calJeffreysDistanceMesh(Gaussian *ge){
  //  cout << "GaussianMixture::calJeffreysDistanceMesh"<<endl;
  double kldiv12 = calKLDivergenceMesh(ge);
  double kldiv21 = ge->calKLDivergenceMesh(this);
  //    cout << "dbg kl1:"<<kldiv12<<" "<<kldiv21<<endl;
  return kldiv12+kldiv21;
}
/////////////////// GaussianMixture /////////////////////

GaussianMixture::GaussianMixture(){
}
int GaussianMixture::pushMixtureElement(double inPi, Gaussian inElem){
  pi.push_back(inPi);
  gauss.push_back(inElem);
  return gauss.size();
}

double GaussianMixture::calProb(ublas::vector<double> x){
  vector<double>::iterator itrPi = pi.begin();
  vector<Gaussian>::iterator itrG = gauss.begin();
  double p = 0;
  for(;itrPi!=pi.end();itrPi++, itrG++){
    p += (*itrPi)*(*itrG).calProb(x);
  }
  return p;
}

double GaussianMixture::calKLDivergenceMesh(GaussianMixture *gm){
  //  cout << "GaussianMixture::calKLDivergenceMesh"<<endl;
  double kld = 0.0;
  for(int i=0;i<(int)mesh.size();i++){
    double p1n = mesh[i];
    double p2n = gm->getMesh(i);
    if(p1n>1e-10 && p2n>1e-10)
      kld += p1n*log(p1n/p2n);
    else
      cerr<< "probability is too small."<<endl;
  }
  return kld;
}
double GaussianMixture::calJeffreysDistanceMesh(GaussianMixture *gm){
  //  cout << "GaussianMixture::calJeffreysDistanceMesh"<<endl;
  double kldiv12 = calKLDivergenceMesh(gm);
  double kldiv21 = gm->calKLDivergenceMesh(this);
  //    cout << "dbg kl1:"<<kldiv12<<" "<<kldiv21<<endl;
  return kldiv12+kldiv21;
}

double GaussianMixture::setProbabilityMesh(ublas::vector<double> init,
					   ublas::vector<double> stepwidth,
					   ublas::vector<int> n_steps,
					   bool xyzToSpherical,
					   double laplace_delta){
  //  cout << "GaussianMixture::setProbabilityMesh"<<endl;
  int dimension = init.size();
  ublas::vector<double> var = init;
  ublas::vector<int> cur_step(dimension);
  for(int i=0; i<dimension; i++) cur_step[i]=0;
  int sum_n_steps=1;
  ublas::vector<int>::iterator itr_nsteps;
  for(itr_nsteps=n_steps.begin();
      itr_nsteps!=n_steps.end(); itr_nsteps++){
    sum_n_steps *= *itr_nsteps;
  }
  if (sum_n_steps == 0) {cout <<"sum_n_steps = 0!!"<<endl; exit(1);}
  mesh = ublas::vector<double>(sum_n_steps);
  for(int i=0;i<sum_n_steps;i++){mesh[i]=0.0;}
  mesh_sum = 0.0;
  for(int i=0; i<sum_n_steps; i++){
    int mlt = 1;
    for(int d=0; d<n_steps.size();d++){
      cur_step[d] = (i/mlt)%n_steps[d];
      var[d] = init[d]+cur_step[d]*stepwidth[d];
      mlt *= n_steps[d];
    }
    if(xyzToSpherical){
      Coord var_c(var[0],var[1],var[2]);
      if(var[0]==0 || var[1]==0 || var[2]==0) continue;
      var[0] = var_c.distTo(Coord(0,0,0));
      var[1] = var_c.calIgTheta();
      var[2] = var_c.calIgPhi();
    }
    double prob = calProb(var)+laplace_delta;
    mesh[i] = prob;
    mesh_sum += prob;
  }
  for(int i=0;i<sum_n_steps;i++){mesh[i]/=mesh_sum;}
  return 0;
}

