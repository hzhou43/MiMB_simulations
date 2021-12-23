//#include <iostream>
//#include <iomanip>
#include <cstdio>
#include "timebase.hpp"

TimeBase::TimeBase(int n, std::string *jobs){
  names=jobs;
  njob=n;
  mstart=new double[njob+1];
  msum=new double[njob+1];
#ifdef _OPENMP
  cstart=new double[njob+1];
  csum=new double[njob+1];
#endif
  start(njob);
}

TimeBase::~TimeBase(){
  delete[] mstart;
  delete[] msum;
#ifdef _OPENMP
  delete [] cstart;
  delete [] csum;
#endif
}

void TimeBase::init(){
  for (int i=0;i<njob+1;i++){
	msum[i]=0.0;
  }
#ifdef _OPENMP
  for (int i=0;i<njob+1;i++){
	csum[i]=0.0;
  }
#endif
}

void TimeBase::start(int i){
  double val;
  val = (double) clock() / (double) CLOCKS_PER_SEC;
  if (i<=njob){
    mstart[i]=val;
#ifdef _OPENMP
    val = omp_get_wtime();
    cstart[i]=val;
#endif
  }
}

void TimeBase::starts(int n, ...){
  va_list arguments;
  va_start ( arguments, n );
  int i;
  for (i = 0; i < n; i++ ){
    start(va_arg ( arguments, int ));
  }
}

void TimeBase::stop(int i){
  double val;
  val = (double) clock() / (double) CLOCKS_PER_SEC;
  if (i<=njob){
    msum[i]+=(val-mstart[i]);
#ifdef _OPENMP
    val = omp_get_wtime();
    csum[i]+=(val-cstart[i]);
#endif
  }
}

void TimeBase::stops(int n, ...){
  va_list arguments;
  va_start ( arguments, n );
  int i;
  for (i = 0; i < n; i++ ){
    stop(va_arg ( arguments, int ));
  }
}

void TimeBase::report(FILE *fp){
  int i;
  stop(njob);
  for (i=0;i<njob;i++){
    //std::cout<<std::setw(16)<<names[i]<<":"<<msum[i]<<msum[i]/msum[njob]<<std::endl;
    fprintf(fp,"#Time %-8s%12.3f (Sec) %8.3f (%%)\n",names[i].c_str(),msum[i],msum[i]/msum[njob]*100);
  }
  fprintf(fp,"#Time %-8s%12.3f (Sec) %8.3f (%%)\n","Total:",msum[i],msum[i]/msum[njob]*100);
#ifdef _OPENMP
  for (i=0;i<njob;i++){
    fprintf(fp,"#TimeOMP  %-8s%12.3f (Sec) %8.3f (\%)\n",names[i].c_str(),csum[i],csum[i]/csum[njob]*100);
    //std::cout<<std::setw(16)<<names[i]<<":"<<csum[i]<<csum[i]/csum[njob]<<std::endl;
  }
  fprintf(fp,"#TimeOMP  %-8s%12.3f (Sec) %8.3f (\%)\n","Total:",csum[i],csum[i]/csum[njob]*100);
#endif
}
