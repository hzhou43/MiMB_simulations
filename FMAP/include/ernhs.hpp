#ifndef ernhs_hpp__
#define ernhs_hpp__ 1

#include "ernbase.hpp"

class ErnHS:public ErnBase {
  public:
    ///Constructor eps(r-sig)^2 for r<d
    ErnHS(double c, double eps, double sig):ErnBase(c){
      //rc2=c*c;
      meps=eps;
      msig=sig;
      msig2=msig*msig;
    }
    ErnHS(double l, double c, double eps, double sig, double offset):ErnBase(l,c){
      //rc2=c*c;
      meps=eps;
      msig=sig;
      msig2=msig*msig;
    }
    ~ErnHS(){}
    
    ///Return ern with distance d^2.
    double ern2(double d2);

    ///Return virial presure with distance d^2.
    double vir2(double d2);

    private:
      double rc2;
      double msig2;
      double meps;
      double msig;
      ///actual energy function
      double hs(double d2);
};
#endif //ernhs_hpp__
