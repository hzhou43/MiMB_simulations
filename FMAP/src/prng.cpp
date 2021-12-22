#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
//#include <gsl/gsl_rng.h>
#include "prng.hpp"

using namespace std;

//PrngBase
int PrngBase::ranI(int up){
  return (int)floor(ran()*up);
}

double PrngBase::check(int n){
  double sum=0.0;
  for (int i=0;i<n;i++){
    sum+=this->ran();
  }
  cout<< "Check random genrator, expected: "<< n*0.5 <<", get: "<< sum <<"."<<endl;
  return sum;
}

//PrngSys
void PrngSys::init(){
  //cout<<"init_ran"<<endl;
  srand(seed);
}

double PrngSys::ran(){
  return (double)rand()/RAND_MAX;
}
/*
//PrngGSL
PrngGSL::~PrngGSL(){
  if (prng){
    gsl_rng_free((gsl_rng *)prng);
    prng=NULL;
  }
}

void PrngGSL::init(){
  //cout<<"gsl init_ran"<<endl;
  prng = (void *)gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set((gsl_rng *)prng,seed);
}

double PrngGSL::ran(){
  //cout<<"gsl ran"<<endl;
  return gsl_rng_uniform((gsl_rng *)prng);
}

int PrngGSL::ranI(int up){
  return (int)gsl_rng_uniform_int((gsl_rng *)prng,up);
}
*/
//PrngNR
struct RANP {
        long id;
        int in;
        int inp;
        long ma[56];
};

typedef struct RANP *RANPARM;

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

RANPARM ran1_init(int seed){
        RANPARM ranp;
        if(!(ranp = (RANPARM)malloc(sizeof(struct RANP)))) {
                return 0;
        }
        ranp->id=seed;
        return ranp;
}

/// return a random number between [0:1.0)
double ran1(RANPARM ranp){
        long mj,mk;
        int i,ii,k;

        if (ranp->id <0 ){
                mj=MSEED-(ranp->id < 0 ? -ranp->id : ranp->id);
                mj %= MBIG;
                ranp->ma[55]=mj;
                mk=1;
                for (i=1;i <= 54;i++) {
                        ii=(21*i) % 55;
                        ranp->ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ranp->ma[ii];
                }
                for (k=1;k <= 4;k++)
                        for (i=1;i <= 55;i++) {
                                ranp->ma[i] -= ranp->ma[1+(i+30) % 55];
                                if (ranp->ma[i] < MZ) ranp->ma[i] += MBIG;
                        }
                ranp->in=0;
                ranp->inp=31;
                ranp->id=1;
        }
        if (++ranp->in == 56) ranp->in=1;
        if (++ranp->inp == 56) ranp->inp=1;
        mj=ranp->ma[ranp->in]-ranp->ma[ranp->inp];
        if (mj < MZ) mj += MBIG;
        ranp->ma[ranp->in]=mj;
        return mj*FAC;
}

void ran1_free(RANPARM *ranp){
        assert(ranp && *ranp);
        if (ranp && *ranp){
                free(*ranp);
                *ranp=NULL;
        }
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

PrngNR::~PrngNR(){
  if (prng){
    ran1_free((RANPARM *)&prng);
  }
}

void PrngNR::init(){
  if (seed==0){
    seed=23410981;
  }
  if (seed>0){
    seed=-seed;
  }
  prng=(void *)ran1_init(seed);
}

double PrngNR::ran(){
  return ran1((RANPARM)prng);
}

// fake PrngGSL with PrngNR
PrngGSL::~PrngGSL(){
  if (prng){
    ran1_free((RANPARM *)&prng);
  }
}

void PrngGSL::init(){
  if (seed==0){
    seed=23410981;
  }
  if (seed>0){
    seed=-seed;
  }
  prng=(void *)ran1_init(seed);
}

double PrngGSL::ran(){
  return ran1((RANPARM)prng);
}
