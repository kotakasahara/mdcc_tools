#include "Kprml.h"

Kprml::Kprml(){
}

int Kprml::setup(int argn, char* argv[]){
  string fnCfg;
  if(argn<2){
    cerr << "Usage: kprml [mode]"<<endl;
    cerr << "------------------------------------"<<endl;
    exit(1);
  }
  cout << "conf.setall\n";
  cfg.setAll(argn,argv);
  if(cfg.fnCfg!=string())
    cfg.setAll(Read(cfg.fnCfg).loadConfig());
  cout <<"/setup\n";
  return 0;
}

int Kprml::mainStream(){
  cout<<"mainstream\n";
  switch(cfg.mode){
  case M_TEST: testMode();   break;
  case M_EM: emMode();       break;
  case M_VBGMM: vbgmmMode(); break;
  case M_VBFULL: vbfullMode(); break;
  case M_CAL_GM_DIST: calGaussianMixtureJeffreysDistanceMode(); break;
  case M_CAL_GE_DIST: calGaussianElementJeffreysDistanceMode(); break;
  default:
    cout <<"mode ga shitei sareteimasen.\n";
    testMode();
    break;
  }
  cout <<"finished.\n";
  return 0;
}

int Kprml::testMode(){
  cout<<"test mode dayo.\n";
  ublas::vector<ublas::vector<double> > dataTable;
  cout << "read table"<<endl;
  if(cfg.formatDataTable == DATA_TABLE){
    cout << "read tsv table"<<endl;
    dataTable=Read(cfg.fnDataTable).loadDataTable(cfg.featureDef, cfg.dataSkip, cfg.headerSkip, cfg.endFrame);
  }else if(cfg.formatDataTable == DATA_KKTRJ)
    cout << "read kktrajtrans"<<endl;
  dataTable=Read(cfg.fnDataTable).loadKKTrajTrans(cfg.featureDef, cfg.dataSkip, cfg.headerSkip, cfg.endFrame);;
  cout << "dump"<<endl;

  for(int n=0; n<(int)dataTable.size(); n++){
    cout << n ;
    for(int d=0; d<(int)dataTable[n].size(); d++){
      cout << " " << dataTable[n][d];
    }
    cout << endl;
  }
  
  return 0;
}

int Kprml::emMode(){
  cout << "mode em"<<endl;
  ublas::vector<ublas::vector<double> > dataTable;
  cout << "read table"<<endl;
  dataTable=Read(cfg.fnDataTable).loadDataTable(cfg.featureDef, cfg.dataSkip);
  cout << "em"<<endl;
  EMAlgorithm em(dataTable, cfg.fnOutGaussian, cfg.nMixedElement,
		 cfg.maxIteration, cfg.thresholdConvergence,
		 cfg.printInterval,
		 cfg.maxIterationKmeans);
  cout << "em.initialize()"<<endl;
  //  em.initialize();
  em.kmeansRepeat(100);
  cout << "em.iteration()"<<endl;
  em.iteration();
  cout << "output"<<endl;
  
  return 0;
}
int Kprml::vbgmmMode(){
  cout << "mode vbgmm"<<endl;
  ublas::vector<ublas::vector<double> > dataTable;
  cout << "read table"<<endl;
  dataTable=Read(cfg.fnDataTable).loadDataTable(cfg.featureDef,  cfg.dataSkip, cfg.headerSkip);
  cout << "vbgmm"<<endl;
  VBgmm vbgmm(dataTable, cfg.fnOutGaussian, cfg.nMixedElement,
	      cfg.maxIteration, cfg.thresholdConvergence,
	      cfg.printInterval,
	      cfg.vbgmm_hp_beta,
	      cfg.maxIterationKmeans);
  cout << "vbgmm.main()"<<endl;
  vbgmm.main();
  return 0;
}

int Kprml::vbfullMode(){
  cout << "mode vbfull"<<endl;
  ublas::vector<ublas::vector<double> > dataTable;
  if(cfg.formatDataTable == DATA_TABLE){
    cout << "Read tsv table"<<endl;
    dataTable=Read(cfg.fnDataTable).loadDataTable(cfg.featureDef, cfg.dataSkip, cfg.headerSkip);
  }else if(cfg.formatDataTable == DATA_KKTRJ){
    cout << "Read mDCC binary data table"<<endl;
    dataTable=Read(cfg.fnDataTable).loadKKTrajTrans(cfg.featureDef, cfg.dataSkip, cfg.headerSkip);
  }

  VBfull vbfull(dataTable, 
		cfg.fnOutGaussian,
		cfg.nMixedElement,
		cfg.maxIteration, cfg.thresholdConvergence,
		cfg.printInterval,
		cfg.vbfull_alpha_0,
		cfg.vbfull_w_0_coef,
		cfg.vbfull_ny_0,
		cfg.maxIterationKmeans);
  //cout << "vbfull.main()"<<endl;
  vbfull.main();
  return 0;
}

int Kprml::calGaussianMixtureJeffreysDistanceMode(){
  cout << "mode calGaussianMixtureJeffreysDistance"<<endl;

  int n_pair;
  int cur;
  double max_dist = 0.0;

  map<string,GaussianMixture> gm;
  int dim = cfg.klGridStepwidth.size();
  ublas::vector<double> v_init(dim);
  ublas::vector<double> v_stepwidth(dim);
  ublas::vector<int> v_nsteps(dim);
  Read(cfg.fnGaussianMixtureInteractions)
    .loadGaussianMixtures(cfg.gaussianMixtureTypeLength,
			  dim, &gm);
  cout << "gm.size:\t"<<gm.size()<<endl;
  for(int i=0;i<dim;i++){
    v_init[i] = cfg.klGridInit[i];
    v_stepwidth[i] = cfg.klGridStepwidth[i];
    v_nsteps[i] = cfg.klGridNSteps[i];
  }
  calGaussianMixtureJeffreysDistance(&gm,
				     v_init, v_stepwidth, v_nsteps,
				     cfg.flgKlGridConvSpherical);

  return 0;
}

int Kprml::calGaussianMixtureJeffreysDistance(map<string,GaussianMixture> *gm,
					      ublas::vector<double> v_init,
					      ublas::vector<double> v_stepwidth,
					      ublas::vector<int> v_nsteps,
					      bool xyzToSpherical){
  map<string,GaussianMixture>::iterator itr_gm_1;
  map<string,GaussianMixture>::iterator itr_gm_2;
  set<string> flg_mesh_cal;
  int n_pair;
  n_pair = (gm->size() * gm->size() - gm->size())/2;
  int cur = -1;
  int cur1 = -1;
  double max_dist=0.0;

  for(itr_gm_1 = gm->begin();
      itr_gm_1 != gm->end(); itr_gm_1++){
    cur1++;
    if(cfg.klCalcPairBegin1>-1 && cur1<cfg.klCalcPairBegin1) continue;
    if(cfg.klCalcPairEnd1>-1   && cur1 >= cfg.klCalcPairEnd1) break;
    int cur2 = -1;
    for(itr_gm_2 = gm->begin();
	itr_gm_2 != itr_gm_1; itr_gm_2++){
      cur++;
      cur2++;
      if(cfg.klCalcPairBegin>-1 && cur<cfg.klCalcPairBegin) continue;
      if(cfg.klCalcPairEnd>-1   && cur >= cfg.klCalcPairEnd ) break;
      if(cfg.klCalcPairBegin2>-1 && cur2<cfg.klCalcPairBegin2) continue;
      if(cfg.klCalcPairEnd2>-1   && cur2 >= cfg.klCalcPairEnd2) continue;
      
      //      if(cur%10000 == 0) cout<<"tri_atom:\t"<<cur<<" / "<<n_pair<<endl;
      
      if(flg_mesh_cal.find(itr_gm_1->first)==flg_mesh_cal.end()){
	itr_gm_1->second.setProbabilityMesh(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
	flg_mesh_cal.insert(itr_gm_1->first);
      }
      if(flg_mesh_cal.find(itr_gm_2->first)==flg_mesh_cal.end()){
	itr_gm_2->second.setProbabilityMesh(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
	flg_mesh_cal.insert(itr_gm_2->first);
      }

      double dist = itr_gm_1->second.calJeffreysDistanceMesh(&itr_gm_2->second);
      if(cfg.maxDistanceForGaussianMixtures < 1e-10 ||
	 dist <= cfg.maxDistanceForGaussianMixtures)
	cout << "dist:\t"<<cur1<<"\t"<<cur2<<"\t"<<itr_gm_1->first<<"\t"<<itr_gm_2->first<<"\t"
	     << dist << endl;
      if(dist>max_dist) max_dist=dist;
    }
  }
  cout << "max dist:\t"<< max_dist << endl;
  return 0;
}

int Kprml::calGaussianElementJeffreysDistanceMode(){
  cout << "mode calGaussianElementJeffreysDistance"<<endl;

  int n_pair;
  int cur;
  double max_dist = 0.0;

  map<string,GaussianMixture> gm;
  int dim = cfg.klGridStepwidth.size();
  ublas::vector<double> v_init(dim);
  ublas::vector<double> v_stepwidth(dim);
  ublas::vector<int> v_nsteps(dim);
  Read(cfg.fnGaussianMixtureInteractions)
    .loadGaussianMixtures(cfg.gaussianMixtureTypeLength,
			  dim, &gm);
  cout << "gm.size:\t"<<gm.size()<<endl;
  for(int i=0;i<dim;i++){
    v_init[i] = cfg.klGridInit[i];
    v_stepwidth[i] = cfg.klGridStepwidth[i];
    v_nsteps[i] = cfg.klGridNSteps[i];
  }
  calGaussianElementJeffreysDistance(&gm,
				     v_init, v_stepwidth, v_nsteps,
				     cfg.flgKlGridConvSpherical);

  return 0;
}

int Kprml::calGaussianElementJeffreysDistance(map<string,GaussianMixture> *gm,
					      ublas::vector<double> v_init,
					      ublas::vector<double> v_stepwidth,
					      ublas::vector<int> v_nsteps,
					      bool xyzToSpherical){
  map<string,GaussianMixture>::iterator itr_gm_1;
  set<int> flg_mesh_cal;
  int n_pair;
  n_pair = (gm->size() * gm->size() - gm->size())/2;
  int cur = -1;
  int cur1 = -1;
  double max_dist=0.0;
  
  for(itr_gm_1 = gm->begin();
      itr_gm_1 != gm->end(); itr_gm_1++){
    vector<Gaussian>::iterator itr_ge_1;
    for(itr_ge_1 = itr_gm_1->second.getGaussBegin();
	itr_ge_1 != itr_gm_1->second.getGaussEnd(); itr_ge_1++){
      cur1++;
      //      cout << "cur1 " << cur1 << ":" << cfg.klCalcPairBegin1 << " : " << cfg.klCalcPairEnd1 <<endl;
      if(cfg.klCalcPairBegin1>-1 && cur1<cfg.klCalcPairBegin1) continue;
      if(cfg.klCalcPairEnd1>-1   && cur1 >= cfg.klCalcPairEnd1) break;
      int cur2 = -1;
      Gaussian ge1;
      Gaussian *p_ge1 = &ge1;
      if(cfg.modeSaveMemory==1){
	cout << "a\n";
	/* save memory */
	ge1 = (*itr_ge_1);
	ge1.setProbabilityMesh(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
	p_ge1 = &ge1;
      }else if(cfg.modeSaveMemory==2){
	if(flg_mesh_cal.find((*itr_ge_1).getId())==flg_mesh_cal.end()){
	  //	  cout << "# cal mesh sum "<<(*itr_ge_1).getId()<<endl;
	  (*itr_ge_1).setProbabilityMeshSum(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
	  p_ge1 = &*itr_ge_1;
	}
      }else{
	cout << "b\n";
	/* save CPU time */
	if(flg_mesh_cal.find((*itr_ge_1).getId())==flg_mesh_cal.end()){
	  (*itr_ge_1).setProbabilityMesh(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
	  flg_mesh_cal.insert((*itr_ge_1).getId());
	  //	cout << "# cal mesh "<<(*itr_ge_1).getId()<<endl;
	  p_ge1 = &*itr_ge_1;
	}
      }

      map<string,GaussianMixture>::iterator itr_gm_2;
      map<string,GaussianMixture>::iterator itr_gm_end = itr_gm_1;
      itr_gm_end++;
      for(itr_gm_2 = gm->begin();
	  itr_gm_2 != itr_gm_end; itr_gm_2++){
	vector<Gaussian>::iterator itr_ge_2;
	vector<Gaussian>::iterator itr_ge_end;
	if(itr_gm_1 == itr_gm_2)
	  itr_ge_end = itr_ge_1;
	else
	  itr_ge_end = itr_gm_2->second.getGaussEnd();
	for(itr_ge_2 = itr_gm_2->second.getGaussBegin();
	    itr_ge_2 != itr_gm_2->second.getGaussEnd(); itr_ge_2++){
	  cur++;
	  cur2++;
	  //	  cout << "cur2 " << cur2 << ":" << cfg.klCalcPairBegin2 << " : " << cfg.klCalcPairEnd2 <<endl;
	  if(cfg.klCalcPairBegin>-1 && cur<cfg.klCalcPairBegin) continue;
	  if(cfg.klCalcPairEnd>-1   && cur >= cfg.klCalcPairEnd ) break;
	  if(cfg.klCalcPairBegin2>-1 && cur2<cfg.klCalcPairBegin2) continue;
	  if(cfg.klCalcPairEnd2>-1   && cur2 >= cfg.klCalcPairEnd2) continue;
      //      if(cur%10000 == 0) cout<<"tri_atom:\t"<<cur<<" / "<<n_pair<<endl;
	  //	  cout<<"tri_atom:\t"<<cur<<" / "<<n_pair<<endl;
	  
	  double maharanobis_dist_12 = (*itr_ge_1).calGaussMaharanobisDist((*itr_ge_2).getMu());
	  double maharanobis_dist_21 = (*itr_ge_2).calGaussMaharanobisDist((*itr_ge_1).getMu());
	  //	  cout << "# mahara d "<<maharanobis_dist_12 << " " << maharanobis_dist_21 <<endl;
	  cout << "m_dist:\t"<<(*itr_ge_1).getId()<<"\t"<<(*itr_ge_2).getId()<<"\t" // <<itr_gm_1->first<<"\t"<<itr_gm_2->first<<"\t"
	       << maharanobis_dist_12 << "\t" << maharanobis_dist_21 << endl;
	  if(cfg.klMaxMaharanobisDist < 1e-10 ||
	     (maharanobis_dist_12 <= cfg.klMaxMaharanobisDist ||
	      maharanobis_dist_21 <= cfg.klMaxMaharanobisDist)){
	    Gaussian ge2;
	    Gaussian *p_ge2 = &ge2;
	    if(cfg.modeSaveMemory==1){
	      /* save memory */
	      ge2 = (*itr_ge_2);
	      ge2.setProbabilityMesh(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
	      p_ge2 = &ge2;
	    }else if(cfg.modeSaveMemory==2){
	      if(flg_mesh_cal.find((*itr_ge_2).getId())==flg_mesh_cal.end()){
		(*itr_ge_2).setProbabilityMeshSum(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
		//	cout << "# cal mesh "<<(*itr_ge_1).getId()<<endl;
		p_ge2 = &*itr_ge_2;
	      }
	    }else{
	      /* save CPU time */
	      if(flg_mesh_cal.find((*itr_ge_2).getId())==flg_mesh_cal.end()){
		(*itr_ge_2).setProbabilityMesh(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta);
		flg_mesh_cal.insert((*itr_ge_2).getId());
		//	cout << "# cal mesh "<<(*itr_ge_1).getId()<<endl;
		p_ge2 = &*itr_ge_2;
	      }
	    }
	    
	    double dist=0;
	    double kl21=0;
	    double kl12=0;
	    if(cfg.modeSaveMemory==2){
	      kl12 = p_ge1->calKlDivergence(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta, p_ge2);
	      kl21 = p_ge2->calKlDivergence(v_init,v_stepwidth,v_nsteps,xyzToSpherical, cfg.klCalcLaplaceDelta, p_ge1);
	      dist = kl12 + kl21;
	    }else{
	      dist = p_ge1->calJeffreysDistanceMesh(p_ge2);
	      //	    cout << "# jeff_dist " <<dist <<endl;
	    }
	    if(cfg.maxDistanceForGaussianMixtures < 1e-10 ||
	       dist <= cfg.maxDistanceForGaussianMixtures)
	      cout << "dist:\t"<<p_ge1->getId()<<"\t"<<p_ge2->getId()<<"\t" // <<itr_gm_1->first<<"\t"<<itr_gm_2->first<<"\t"
		   << dist << "\t" << kl12 << "\t" <<  kl21 << endl;
	    if(dist>max_dist) max_dist=dist;
	    if(cfg.modeSaveMemory>=1){
	      (*itr_ge_2).clearProbabilityMesh();
	      set<int>::iterator itrTmp = flg_mesh_cal.find((*itr_ge_2).getId());
	      if(itrTmp!=flg_mesh_cal.end()) flg_mesh_cal.erase(itrTmp);
	    }
	  }
	}
      }
      if(cfg.modeSaveMemory>=1){
	(*itr_ge_1).clearProbabilityMesh();
	set<int>::iterator itrTmp = flg_mesh_cal.find((*itr_ge_1).getId());
	if(itrTmp!=flg_mesh_cal.end()) flg_mesh_cal.erase(itrTmp);
      }
    }
  }
  cout << "max dist:\t"<< max_dist << endl;
  return 0;
}
