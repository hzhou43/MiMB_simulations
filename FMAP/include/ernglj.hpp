#ifndef ernglj_hpp__
#define ernglj_hpp__ 1

#include "ernbase.hpp"

class ErnGLJ:public ErnBase {
  public:
    ///Constructor 4e(s^12/(r-o)^12-s^6/(r-o)^6)
    ErnGLJ(double c, double eps, double sig, double offset):ErnBase(c){
      rc2=c*c;
      meps=eps;
      msig=sig;
      moffset=offset;
    }
    ErnGLJ(double l, double c, double eps, double sig, double offset):ErnBase(l,c){
      rc2=c*c;
      meps=eps;
      msig=sig;
      moffset=offset;
    }
    ~ErnGLJ(){}
    
    ///Return ern with distance d^2.
    double ern2(double d2);

    ///Return virial presure with distance d^2.
    double vir2(double d2);

    private:
      double rc2;
      double meps;
      double msig;
      double moffset;
      ///actual energy function
      double glj(double d2);
};
#endif //ernglj_hpp__
