#include <cmath>
#include "ernhs.hpp"

double ErnHS::hs(double d2){
  if (d2<msig2){
    //double r=sqrt(d2);
    //double dr=r-msig;
    //return meps*dr*dr;
    return meps;
  }else{
    return 0.0;
  }
}

double ErnHS::ern2(double d2){
  return hs(d2);
}

//Fake useless
double ErnHS::vir2(double d2){
    return 0.0;
}
