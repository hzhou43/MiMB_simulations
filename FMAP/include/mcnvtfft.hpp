#ifndef mcnvtfft_hpp__
#define mcnvtfft_hpp__ 1

#include <assert.h>
#include <iostream>
#include "mcnvtbase.hpp"
#include "fftw.h"

///Base class for MC nvt.
class MCNVTFFT:public MCNVTBase {
  public:
    MCNVTFFT(int n, double l, double t, bool r, PrngBase *p, int c, ErnBase *e, int nf, double vs);
    ~MCNVTFFT();
    
    void buildR(PrngBase *p);
    void calc();
    double getval(int i, int j, int k);
  protected:
    int nfold;
    int ngrid;
    double spacing;
    double vscl;
  
  private:
    fftw_real *grdrec,*grdlig;
    double *xf,*yf,*zf,*ligx,*ligy,*ligz;
    double *coef, *coef2;
    int *xi,*yi,*zi;
};
#endif //mcnvtfft_hpp__
