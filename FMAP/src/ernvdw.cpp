#include "ernvdw.hpp"
#include "ljs.h"

double ErnVdwLfs::ern2(double d2){
  if (d2<rc2){
    //return LjLfsSlow(d2,rc2);
    return LjShfSlow(d2,rc2);
  }else{
    return 0.0;
  }
}

double ErnVdwLfs::vir2(double d2){
  if (d2<rc2){
    //return LjLfsWRc(d2,rc2);
    return LjShfWRc(d2,rc2);
  }else{
    return 0.0;
  }
}

double ErnVdwShf::ern2(double d2){
  if (d2<rc2){
    return LjShfSlow(d2,rc2);
  }else{
    return 0.0;
  }
}

double ErnVdwShf::vir2(double d2){
  if (d2<rc2){
    return LjShfWRc(d2,rc2);
  }else{
    return 0.0;
  }
}

