#include <iostream>
#include "define.h"
#include "Kprml.h"
using namespace std; 

int main(int argn, char* argv[]){
  cout << ABOUT_ME <<endl;
  Kprml kprml;
  //  cout << "ababa "<<argn<<endl;
  kprml.setup(argn,argv);
  kprml.mainStream();
  return 0;
}
