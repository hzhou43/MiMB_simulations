#ifndef ernvdw_hpp__
#define ernvdw_hpp__ 1

#include "ernbase.hpp"

class ErnVdwLfs:public ErnBase {
  public:
    ///Constructor
    ErnVdwLfs(double c):ErnBase(c){
      rc2=c*c;
    }
    ErnVdwLfs(double l, double c):ErnBase(l,c){
      rc2=c*c;
    }
    ~ErnVdwLfs(){}
    
    ///Return ern with distance d^2.
    double ern2(double d2);

    ///Return virial presure with distance d^2.
    double vir2(double d2);

    private:
      double rc2;
};

class ErnVdwShf:public ErnBase {
  public:
    ///Constructor
    ErnVdwShf(double c):ErnBase(c){
      rc2=c*c;
    }
    ErnVdwShf(double l, double c):ErnBase(l,c){
      rc2=c*c;
    }
    ~ErnVdwShf(){}

    ///Return ern with distance d^2.
    double ern2(double d2);

    ///Return virial presure with distance d^2.
    double vir2(double d2);

    private:
      double rc2;
};

#endif //ernvdw_hpp__
