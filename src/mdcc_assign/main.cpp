#include <iostream>
#include "define.h"
#include "PliGauss.h"
using namespace std; 

int main(int argn, char* argv[]){
  clock_t timestart;
  timestart=clock();
  cout << ABOUT_ME <<endl;
  PliGauss pg;
  //  cout << "ababa "<<argn<<endl;
  pg.setup(argn,argv);
  pg.mainStream();
  cout<<"time "<<(clock()-timestart)/CLOCKS_PER_SEC<<endl;
  return 0;
}
