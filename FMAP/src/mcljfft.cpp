#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "mcnvtfft.hpp"
#include "ernvdw.hpp"
#include "timebase.hpp"
#include "barhist.h"
#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif
int mcljboe(int argc, char** argv);
double mcljrun(int N, double L, double T, double rc, double dr, int nCycles, int nEq, int fs, int ts, int Seed);

#ifdef __EMSCRIPTEN__
#ifdef __cplusplus
extern "C" {
#endif

EM_JS(void, call_alert, (int x), {
  //alert('hello world!');
  console.log('Progress: ' + x + ' %');
  //var t=document.getElementById('sim_prog');
  //console.log('old:'+t.value);
  //t.value=x;
  //console.log('new:'+t.value);
  //throw 'all done';
});


EMSCRIPTEN_KEEPALIVE double mcLJ(int N, double L, double T, double rc, int nCycles, int nEq, int fs, int ts, int Seed){
	double dr=0.3;
	return mcljrun(N,L,T,rc,dr,nCycles,nEq,fs,ts,Seed);	
}

#ifdef __cplusplus
}
#endif
#endif

int main(int argc, char** argv){
	return mcljboe(argc, argv);
}

int mcljboe(int argc, char** argv){
  // Here we parse the command line arguments
  int N=216;
  double L=8,T=1.0,rc=2.5, dr=0.1;
  int nCycles = 10, nEq=1000;
  int fs=1,ts=1;
  int Seed = 23410981;
  double vscl=1.06;
  int scl=10;
 
  int i,j;
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-L")) L=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ne")) nEq = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fs = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ts")) ts = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-vscl")) vscl = atof(argv[++i]);
    else if (!strcmp(argv[i],"-scl")) scl = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-seed"))
      Seed = atoi(argv[++i]);
    else {
      fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",
                 argv[i]);
      exit(-1);
    }
  }
  mcljrun(N,L,T,rc,dr,nCycles,nEq,fs,ts,Seed);
  return 0; 
}

double mcljrun(int N, double L, double T, double rc, double dr, int nCycles, int nEq, int fs, int ts, int Seed){
  int i,j;
  double vscl=1.06;
  int scl=10;
  std::string jobs[]={"Insert:"};
  TimeBase *tp=new TimeBase(1,jobs);
  PrngBase *pb=new PrngGSL(Seed);
  PrngBase *pb2=new PrngGSL(Seed+1);
  PrngBase *pd=new PrngGSL(Seed);
  PrngBase *pd2=new PrngGSL(Seed+1);
  ErnBase  *eb=new ErnVdwLfs(rc); 
  MCNVTFFT *nvt=new MCNVTFFT(N,L,T,true,pb,N,eb,scl,vscl);
  tp->init();
  nvt->init_cube();
  nvt->setdr(dr);
  nvt->acceptratio(true);
  nvt->stat(false);
  nCycles+=nEq;
  double e_sum=0.0,w_sum=0.0;
  long nc=0;
  
  int ng=(int)(L)*scl;
  nvt->buildR(pb2);

  double beta=nvt->beta();
  double we;
  
  for (i=1;i<=nCycles;i++){
    nvt->step();
    if (i>nEq && (!(i%fs))){
      e_sum+=nvt->ern(false);
      tp->start(0);
      for (j=0;j<ts;j++){
	nvt->calc();
	int gi,gj,gk;
	for (gi=0;gi<ng;gi++){
	  for (gj=0;gj<ng;gj++){
	    for (gk=0;gk<ng;gk++){
	      we=nvt->getval(gi,gj,gk);
	      w_sum+=exp(-we*beta);
	    }
	  }
	}
        //we=nvt->widom(pb2);
      }
      tp->stop(0);
      nc++;
    }
#ifdef __EMSCRIPTEN__
    if(!(i%fs)){
	call_alert(i*100/nCycles);
    }
#endif
  }
  double exchem=-log(w_sum/nc/ts/ng/ng/ng);
  fprintf(stdout,"#Density: %.8lf acc: %.8lf avE: %.8lf wE: %.8lf\n",nvt->numdensity(), nvt->acceptratio(false), e_sum/nc/N*beta,exchem);
  fflush(stdout);
  tp->report(stdout);
  delete pb;
  delete pb2;
  delete eb;
  delete nvt;
  delete tp;
  //return(-log(w_sum/nc/ts)/beta);
  return exchem;
}
