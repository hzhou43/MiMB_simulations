#ifndef _fftw_h__
#define _fftw_h__ 1

#include <fftw3.h>
#define fftw_real double

#ifdef __cplusplus
extern "C" {
#endif

void fft3d_r2c(int L, int M, int N, fftw_real *a);
void fft3d_c2r(int L, int M, int N, fftw_real *c);
void fft3d_add_inplace(int L, int M, int N,fftw_real *a, fftw_real *b);

///helper
void setzero(int L, int M, int N, fftw_real *a);
double blindsum(int L, int M, int N, fftw_real *a);
#ifdef __cplusplus
}
#endif

#endif  /* _fftw_h__ */
