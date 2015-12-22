#include "Clustering.h"

Cluster::Cluster(int inId){
  id=inId;
}
Cluster::Cluster(int inId, int inCmpt){
  id=inId;
  insert(inCmpt);
}
Cluster::Cluster(int inId, set<int> inCmpt){
  id=inId;
  component=inCmpt;
}
bool Cluster::insert(int inC){
  return component.insert(inC).second;
}
int Cluster::merge(Cluster c){
  set<int> cCmpt=c.getComponent();
  set<int>::iterator itr;
  for(itr=cCmpt.begin(); itr!=cCmpt.end(); itr++){
    component.insert(*itr);
  }
  return component.size();
}

Clustering::Clustering(){
}
Clustering::Clustering(map<pair<int,int>,double> inCmptDist){
  cmptDist=inCmptDist;
}

int Clustering::insertNewCluster(int inCmpt){
  int newId=clusters.size();
  clusters.insert(make_pair(newId, Cluster(newId, inCmpt)));
  clustLink.insert(make_pair(newId,set<int>()));
  return newId;
}

int Clustering::insertComponentDist(pair<int,int> inCmptId, double inDist){
  cmptDist.insert(make_pair(inCmptId,inDist));
  return cmptDist.size();
}

int Clustering::completeLinkage(double maxDist){
  int nlink=initializeClustLinks(maxDist);
  cout<<"nlink="<<nlink<<endl;
  int newClustId=clusters.size();
  while(!clustPair.empty()){
    int c1Id=clustPairDist.begin()->second.first;
    int c2Id=clustPairDist.begin()->second.second;
    while(clustPair.find(make_pair(c1Id,c2Id))==clustPair.end()){
      clustPairDist.erase(clustPairDist.begin());
      c1Id=clustPairDist.begin()->second.first;
      c2Id=clustPairDist.begin()->second.second;
    }
    clustPairDist.erase(clustPairDist.begin());
    clustPair.erase(clustPair.find(make_pair(c1Id,c2Id)));
    Cluster cNew(newClustId);
    map<int,Cluster>::iterator itrC1=clusters.find(c1Id);
    map<int,Cluster>::iterator itrC2=clusters.find(c2Id);
    cNew.merge(itrC1->second);
    cNew.merge(itrC2->second);
    clusters.insert(make_pair(newClustId,cNew));
    clusters.erase(itrC1);
    clusters.erase(itrC2);
    clustLink.insert(make_pair(newClustId,set<int>()));

    map<int,set<int> >::iterator itrC1Link=clustLink.find(c1Id);
    map<int,set<int> >::iterator itrC2Link=clustLink.find(c2Id);
    while(!itrC1Link->second.empty() &&
	  !itrC2Link->second.empty()){
      int c1pId=*(itrC1Link->second.begin());
      int c2pId=*(itrC2Link->second.begin());
      if(c1pId<c2pId){
	removeOldClustInfo(c1Id,c1pId);
      }else if(c1pId>c2pId){
	removeOldClustInfo(c2Id,c2pId);
      }else{
	double newDist=calDistClustPair(cNew,clusters.find(c1pId)->second);
	pair<int,int> newIdPair(c1pId,newClustId);
	clustPairDist.insert(make_pair(newDist,newIdPair));
	clustPair.insert(newIdPair);
	map<int,set<int> >::iterator itrCNewLink=clustLink.find(newClustId);
	itrCNewLink->second.insert(c1pId);
	map<int,set<int> >::iterator itrC1pLink=clustLink.find(c1pId);
	itrC1pLink->second.insert(newClustId);
	removeOldClustInfo(c1Id,c1pId);
	removeOldClustInfo(c2Id,c2pId);
      }
    }
    while(!itrC1Link->second.empty()){
      int c1pId=*(itrC1Link->second.begin());
      removeOldClustInfo(c1Id,c1pId);
    }
    while(!itrC2Link->second.empty()){
      int c2pId=*(itrC2Link->second.begin());
      removeOldClustInfo(c2Id,c2pId);
    }
    clustLink.erase(itrC1Link);
    clustLink.erase(itrC2Link);
    newClustId++;
  }
  return 0;
}

int Clustering::initializeClustLinks(double maxDist){
  map<int,Cluster>::iterator itrC1;
  map<int,Cluster>::iterator itrC2;
  int nlink=0;
  for(itrC1=clusters.begin(); itrC1!=clusters.end();itrC1++){
    for(itrC2=clusters.begin(); itrC2!=itrC1; itrC2++){
      double dist=calDistClustPair(itrC1->second, itrC2->second);
      if(dist>maxDist || dist==-1) continue;
      clustLink.find(itrC1->first)->second.insert(itrC2->first);
      clustLink.find(itrC2->first)->second.insert(itrC1->first);
      nlink++;
      pair<int,int> idpair(itrC2->first,itrC1->first);
      clustPair.insert(idpair);
      clustPairDist.insert(make_pair(dist,idpair));
    }
  }
  return nlink;
}

double Clustering::calDistClustPair(Cluster c1, Cluster c2){
  set<int> cmpt1=c1.getComponent();
  set<int> cmpt2=c2.getComponent();
  double dist=0;
  set<int>::iterator itr1;
  set<int>::iterator itr2;
  for(itr1=cmpt1.begin();itr1!=cmpt1.end();itr1++){
    for(itr2=cmpt2.begin();itr2!=cmpt2.end();itr2++){
      map<pair<int,int>,double>::iterator itrDist;
      pair<int,int> idpair;
      if(*itr1<*itr2) idpair=make_pair(*itr1,*itr2);
      else idpair=make_pair(*itr2,*itr1);
      itrDist=cmptDist.find(idpair);
      if(itrDist==cmptDist.end()){
	//	cout<<"dist n/a "<<idpair.first<<" "<<idpair.second<<endl;
	return -1;
      }
      if(itrDist->second>dist) dist=itrDist->second;
    }
  }
  return dist;
}

int Clustering::removeOldClustInfo(int cId,int cpId){
  map<int,set<int> >::iterator itrCLink=clustLink.find(cId);
  set<int>::iterator itrCpInCLink=itrCLink->second.find(cpId);
  map<int,set<int> >::iterator itrCpLink=clustLink.find(cpId);
  set<int>::iterator itrCInCpLink=itrCpLink->second.find(cId);
  itrCLink->second.erase(itrCpInCLink);
  itrCpLink->second.erase(itrCInCpLink);
  pair<int,int> idpair;
  if(cId<cpId) idpair=make_pair(cId,cpId);
  else idpair=make_pair(cpId,cId);
  set<pair<int,int> >::iterator itrPair=clustPair.find(idpair);
  if(itrPair!=clustPair.end())
    clustPair.erase(itrPair);
  return 0;
}
