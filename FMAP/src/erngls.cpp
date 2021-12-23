#include <cmath>
#include "ernglj.hpp"
#include "ljs.h"

double ErnGLJ::glj(double d2){
  double r=sqrt(d2)-moffset;
  double r2=r*r;
  return LJes(meps,msig,r2);
}

double ErnGLJ::ern2(double d2){
  if (d2<rc2){
    return glj(d2);
  }else{
    return 0.0;
  }
}

//Fake useless
double ErnGLJ::vir2(double d2){
    return 0.0;
}
