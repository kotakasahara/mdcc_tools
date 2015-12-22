#ifndef __CLUSTERING_H__
#define __CLUSTERING_H__

#include <iostream>
#include <vector>
#include <map>
#include <set>
using namespace std;

class Cluster{
 private:
  int id;
  set<int> component;  
 public:
  Cluster(int inId);
  Cluster(int inId, int inCmpt);
  Cluster(int inId, set<int> inCmpt);
  bool insert(int inC);
  int getId(){return id;};
  set<int> getComponent(){return component;};
  int merge(Cluster c);
};

class Clustering{
 private:
  map<pair<int,int>, double> cmptDist;
  map<int,Cluster> clusters;
  map<int,set<int> > clustLink;
  set<pair<int,int> > clustPair;
  multimap<double, pair<int,int> > clustPairDist;
 public:
  Clustering();
  Clustering(map<pair<int,int>,double> inCmptDist);
  int getClustersSize(){return clusters.size();};
  int insertNewCluster(int inCmpt);
  int insertComponentDist(pair<int,int> inCmptId, double inDist);
  int completeLinkage(double maxDist);
  int initializeClustLinks(double maxDist);
  double calDistClustPair(Cluster c1, Cluster c2);
  int removeOldClustInfo(int cId,int cpId);
  map<int,Cluster>::iterator getClustersBegin(){return clusters.begin();};
  map<int,Cluster>::iterator getClustersEnd(){return clusters.end();};
  
};



#endif
