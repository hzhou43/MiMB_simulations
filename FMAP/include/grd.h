#ifndef _grd_h__
#define _grd_h__ 1

#include "fftw.h"

#ifdef __cplusplus
extern "C" {
#endif

int PBC(const int x, const int mx);
void floorfrac(int n, double x[], int xi[], double xf[]);
void val2grd(const int n, const double v[], const int px[], const int py[], const int pz[], const int l, fftw_real* grid);
void pbc1(const int n, const int x[], int px[], const int lc, const int shift);
void grdR(const int l, fftw_real* grid, \
         const int n, int xi[], int yi[], int zi[], double xf[], double yf[], double zf[], \
	 double coef1[],double coef2[],\
         const double rup, const double rlow, const double dx, double (*rdf)(double,double));
void grdL(const int l, fftw_real* grid,\
             const int n, const int xi[], const int yi[], const int zi[] ,\
	     const double xf[], const double yf[], const double zf[], const double coef1[]);
void grdL0(const int l, fftw_real* grid,\
             const int n, const int xi[], const int yi[], const int zi[] ,\
	     const double xf[], const double yf[], const double zf[], const double coef1[]);

#ifdef __cplusplus
}
#endif

#endif  /* _grd_h__ */
