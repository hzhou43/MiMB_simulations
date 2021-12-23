#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <iostream>
#include "mcnvtfft.hpp"
#include "ljs.h"
#include "grd.h"

#define Rcut 3.0
double rdfLfs(double r2,double scl){
  double rc2=Rcut*Rcut;
  //return LjesLfsSlow(scl,6.0/(5.0+scl),r2,rc2);
  return LjesShfSlow(scl,6.0/(5.0+scl),r2,rc2);
}

MCNVTFFT::MCNVTFFT(int n, double l, double t, bool r, PrngBase *p, int c, ErnBase *e, int nf, double vs):MCNVTBase(n,l,t,r,p,c,e),nfold(nf),vscl(vs){
  spacing=1.0/nfold;
  ngrid=length*nfold;
  grdrec=fftw_alloc_real(sizeof(fftw_real)*ngrid*ngrid*(2*(ngrid/2+1)));
  grdlig=fftw_alloc_real(sizeof(fftw_real)*ngrid*ngrid*(2*(ngrid/2+1)));
  coef=new double[cap];
  coef2=new double[cap];
  xi=new int[cap];
  yi=new int[cap];
  zi=new int[cap];
  xf=new double[cap];
  yf=new double[cap];
  zf=new double[cap];
  ligx=new double[cap];
  ligy=new double[cap];
  ligz=new double[cap];
  int i;
  for (i=0;i<cap;i++){
    coef[i]=1.0;
    coef2[i]=vscl;
  }
}

MCNVTFFT::~MCNVTFFT(){
  fftw_free(grdrec);
  fftw_free(grdlig);
  delete[] coef;
  delete[] coef2;
  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] xf;
  delete[] yf;
  delete[] zf;
  delete[] ligx;
  delete[] ligy;
  delete[] ligz;
}

void MCNVTFFT::buildR(PrngBase *p){
  int l=ngrid;
  double x1[1],y1[1], z1[1];
  if (p==NULL){
    x1[0]=0.0;
    y1[0]=0.0;
    z1[0]=0.0;
  }else{
    x1[0]=p->ran();
    y1[0]=p->ran();
    z1[0]=p->ran();
  }
  int x1i[1],y1i[1],z1i[1];
  double x1f[1],y1f[1],z1f[1];
  floorfrac(1,x1,x1i,x1f);
  floorfrac(1,y1,y1i,y1f);
  floorfrac(1,z1,z1i,z1f);
  setzero(l,l,l,grdrec);
  grdR(l,grdrec,1,x1i,y1i,z1i,x1f,y1f,z1f,coef,coef2,Rcut,0.1,spacing,rdfLfs);
  fft3d_r2c(l,l,l,grdrec);
}

void MCNVTFFT::calc(){
  int l=ngrid;
  int i;
  int yoffset=cap, zoffset=cap*2;
  for (i=0;i<num;i++){
    ligx[i]=cords[i]/spacing;
    ligy[i]=cords[yoffset+i]/spacing;
    ligz[i]=cords[zoffset+i]/spacing;
  }
  floorfrac(num,ligx,xi,xf);
  floorfrac(num,ligy,yi,yf);
  floorfrac(num,ligz,zi,zf);
  setzero(l,l,l,grdlig);
  grdL(l,grdlig,num,xi,yi,zi,xf,yf,zf,coef);
  fft3d_r2c(l,l,l,grdlig);
  fft3d_add_inplace(l,l,l,grdrec,grdlig);
  fft3d_c2r(l,l,l,grdlig);
}

double MCNVTFFT::getval(int i, int j, int k){
    int kk=2*(ngrid/2+1);
    return grdlig[kk*(i*ngrid+j)+k];
}
