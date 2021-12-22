#include <cstdlib>
#include <cstdio>
#include "atomsbase.hpp"

/*
double AtomsBase::pbcdst(double d, double boxlen){
    double mi,frac;
    frac=modf(fabs(d)/boxlen,&mi);
    frac=(frac<0.5)?frac:(1.0-frac);
    return frac*boxlen;
}
*/

inline double AtomsBase::pbcdst(double d, double boxlen, double hl){
  if (d>hl){
    d -= boxlen;
  }else if(d<-hl){
    d += boxlen;
  }
  return d;
}

double AtomsBase::distance(int i, int j){
  return sqrt(distance2(i,j));
}

double AtomsBase::distance2(int i, int j){
  double d2=0.0,diff;
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    diff=cords[shift+i]-cords[shift+j];
    d2+=diff*diff;
  }
  return d2;
}

double AtomsBase::distance(int i, int j, double boxlen){
  return sqrt(distance2(i,j,boxlen));
}

double AtomsBase::distance2(int i, int j, double boxlen){
  double d2=0.0,diff;
  int d,shift;
  double hl=boxlen/2.0;
  for (d=0;d<dim;d++){
    shift=d*cap;
    diff=cords[shift+i]-cords[shift+j];
    diff=pbcdst(diff,boxlen,hl);
    d2+=diff*diff;
  }
  return d2;
}

double AtomsBase::distance(int i, const std::valarray<double> &probe){
  return sqrt(distance2(i,probe));
}

double AtomsBase::distance2(int i, const std::valarray<double> &probe){
  double d2=0.0,diff;
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    diff=cords[shift+i]-probe[d];
    d2+=diff*diff;
  }
  return d2;
}

double AtomsBase::distance(int i, const std::valarray<double> &probe, double boxlen){
  return sqrt(distance2(i,probe,boxlen));
}

double AtomsBase::distance2(int i, const std::valarray<double> &probe, double boxlen){
  double d2=0.0,diff;
  int d,shift;
  double hl=boxlen/2.0;
  for (d=0;d<dim;d++){
    shift=d*cap;
    diff=cords[shift+i]-probe[d];
    diff=pbcdst(diff,boxlen,hl);
    d2+=diff*diff;
  }
  return d2;
}

void AtomsBase::showcords(bool txt){
  showcords(stdout,txt);
}

void AtomsBase::showcords(FILE *fp, bool txt){
  int i,d;
  int cnt=0;
  if (txt){
    fprintf(fp,"%d\n",cap);
  }
  for (i=0;i<cap;i++){
    cnt++;
    if (!txt){
      fprintf(fp,"ATOM  %5d  Na  Na   %5d   ",cnt,cnt);
    }
    for (d=0;d<dim;d++){
      if (txt){
	fprintf(fp,"%16.8lf",cords[d*cap+i]);
      }else{
	fprintf(fp,"%8.3lf",cords[d*cap+i]);
      }
    }
    fprintf(fp,"\n");
  }
}

void AtomsBase::save(int i){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    dump[d]=cords[shift+i];
  }
}

void AtomsBase::retrieve(int i){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    cords[shift+i]=dump[d];
  }
}


void AtomsBase::set(int i, const std::valarray<double> &x){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    cords[shift+i]=x[d];
  }
}

void AtomsBase::displace(int i, const std::valarray<double> &x){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    cords[shift+i]+=x[d];
  }
}

void AtomsBase::displace(int i, double x[]){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    cords[shift+i]+=x[d];
  }
}

void AtomsBase::displace(int i, const std::valarray<double> &x, double boxlen){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    cords[shift+i]+=x[d];
    if (cords[shift+i]<0.0) cords[shift+i]+=boxlen;
    if (cords[shift+i]>boxlen) cords[shift+i]-=boxlen;
  }
}

void AtomsBase::displace(int i, double x[], double boxlen){
  int d,shift;
  for (d=0;d<dim;d++){
    shift=d*cap;
    cords[shift+i]+=x[d];
    if (cords[shift+i]<0.0) cords[shift+i]+=boxlen;
    if (cords[shift+i]>boxlen) cords[shift+i]-=boxlen;
  }
}
