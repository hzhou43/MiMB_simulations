#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <iostream>
#include "mcnvtbase.hpp"

void MCNVTBase::step(){
  /// Randomly select a particle
  int i=prng->ranI(num);
  /// calculate displacement
  double erni_old=ern(i);
  double dx[3];
  dx[0] = mdr*(0.5-prng->ran());
  dx[1] = mdr*(0.5-prng->ran());
  dx[2] = mdr*(0.5-prng->ran());
  save(i);
  displace(i,dx,length);
  //std::cout << "dx: " << dx[0] <<" "<< dx[1] <<" "<< dx[2] <<std::endl;
  double erni=ern(i);
  double erndf=(erni-erni_old);
  double pr=prng->ran();
  double ee=exp(-mbeta*erndf);
  if (pr<ee){
    mern+=erndf;
    nacc++;
  }else{
    retrieve(i);
  }
  //std::cout <<"erni_old: " <<erni_old <<" erni: " << erni <<std::endl;
  //std::cout <<"pr: " <<pr <<"ee: " <<ee <<"nacc: " <<nacc <<std::endl;
  ntot++;
}

void MCNVTBase::stat(bool verb){
  MCBoxBase::stat(verb);
  //std::cout<<std::scientific;
  //std::cout<<"Ern: " <<std::setw(16)<<ern() <<"\tVir: " <<std::setw(16)<<vir() <<std::endl;
  mern=ern();
  mvir=vir();
  if (verb){
    printf("%-8s%12.4E\t%-8s%12.4E\n","Ern:",mern,"Vir:",mvir);
  }
}

double MCNVTBase::ern(int i, int j){
  double d2=distance2(i,j,length);
  return ff->ern2(d2);
}

double MCNVTBase::ern(int i){
  int j;
  double e=0.0;
  for (j=0;j<num;j++){
    if (j!=i){
      e+=ern(i,j);
    }
  }
  return e;
}

double MCNVTBase::ern(){
  int i,j;
  double e=0.0;
  for (i=0;i<num-1;i++){
    for (j=i+1;j<num;j++){
	e+=ern(i,j);
    }
  }
  return e;
}

double MCNVTBase::ern(bool redo){
  if (redo){
    return ern();
  }else{
    return mern;
  }
}

double MCNVTBase::ern(const std::valarray<double> &probe){
  int i;
  double e=0.0;
  for (i=0;i<num;i++){
    double d2=distance2(i,probe,length);
    e+=ff->ern2(d2);
  }
  return e;
}

double MCNVTBase::vir(int i, int j){
  double d2=distance2(i,j,length);
  return ff->vir2(d2);
}

double MCNVTBase::vir(int i){
  int j;
  double e=0.0;
  for (j=0;j<num;j++){
    if (j!=i){
      e+=vir(i,j);
    }
  }
  return e;
}

double MCNVTBase::vir(){
  int i,j;
  double e=0.0;
  for (i=0;i<num-1;i++){
    for (j=i+1;j<num;j++){
        e+=vir(i,j);
    }
  }
  return e;
}

double MCNVTBase::vir(const std::valarray<double> &probe){
  int i;
  double e=0.0;
  for (i=0;i<num;i++){
    double d2=distance2(i,probe,length);
    e+=ff->vir2(d2);
  }
  return e;
}

double MCNVTBase::widom(PrngBase *p, ErnBase *f){
  int i;
  double e=0.0;
  if (p==NULL) p=prng;
  if (f==NULL) f=ff;
  dump[0]=(p->ran()-0.5)*length;
  dump[1]=(p->ran()-0.5)*length;
  dump[2]=(p->ran()-0.5)*length;
  for (i=0;i<num;i++){
    double d2=distance2(i,dump,length);
    e+=f->ern2(d2);
  }
  return e;
}

double MCNVTBase::widom(PrngBase *p){
  return widom(p,NULL);
}

double MCNVTBase::widom(){
  return widom(NULL,NULL);
}

double MCNVTBase::widomsbltz(PrngBase *p, int n){
  int i;
  double e=0.0;
  for (i=0;i<n;i++){
    e+=exp(-mbeta*widom(p));
  }
  return e;
}

double MCNVTBase::widomsbltz(int n){
  return widomsbltz(NULL,n);
}

double MCNVTBase::pickdel(PrngBase *p){
  if (p==NULL) p=prng;
  int i=p->ranI(num);
  return ern(i);
}

double MCNVTBase::pickdel(){
  return pickdel(NULL);
}

double MCNVTBase::pickdelsbltz(PrngBase *p,int n){
  int i;
  double e=0.0;
  for (i=0;i<n;i++){
    e+=exp(mbeta*pickdel(p));
  }
  return e;
}

double MCNVTBase::pickdelsbltz(int n){
  return pickdelsbltz(NULL,n);
}
