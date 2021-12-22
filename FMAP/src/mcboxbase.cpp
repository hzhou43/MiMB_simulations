#include <iostream>
#include <cstdlib>
#include "mcboxbase.hpp"

using namespace std;

void MCBoxBase::steps(int n){
  int i;
  for (i=0;i<n;i++){
    step();
  }
}

void MCBoxBase::step(){
  if (prng->ran()>0.5){
    num++;
    nacc++;
  }else{
    num--;
  }
  ntot++;
}

void MCBoxBase::stat(bool verb){
  mbeta=beta();
  if (verb){
    cout << "Volume: " <<vol();
    cout << "\tBeta: "  <<beta();
    cout << "\tDensity: " <<numdensity() <<endl;
    cout << "Accept ratio: " <<acceptratio(false) <<endl;
  }
}

double MCBoxBase::acceptratio(bool init){
  if (init){
    nacc=0;
    ntot=0;
  }
  if (ntot==0){
    return 0.0;
  }else{
    return double(nacc)/(double)(ntot);
  }
}
