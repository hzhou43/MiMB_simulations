/**
 *      @file    bar.c
 *      @author  Sanbo Qin
 *      @brief   Estimate free energy by BAR, OS and EXP.
 */
#include<math.h>
#include "bar.h"

///Generic function, scl(0.5 for 2 stages,1 for 1), drt(forward:1.0,backward:-1.0)
double fIns(double mu, double beta, double parC, double scl, double drt, double (*w)(double,double,double)){
	return w(mu,beta,parC)*exp(-drt*beta*mu*scl);
}

///BAR \f$ w(u)=1/cosh[\beta(\mu-C)/2] \f$.
double wghtBar(double mu, double beta, double parC){
	return 1.0/cosh(beta*((mu-parC)*0.5));
}

///faked, return 1 always.
double wghtOne(double mu, double beta, double parC){
	return 1.0;
}

///BAR insertion
double BarIns(double mu, double beta, double parC){
	return fIns(mu,beta,parC,0.5,1.0,wghtBar);
}

///BAR deletion
double BarDel(double mu, double beta, double parC){
        return fIns(mu,beta,parC,0.5,-1.0,wghtBar);
}

///OS insertion
double OsIns(double mu, double beta){
	return fIns(mu,beta,0.0,0.5,1.0,wghtOne);
}

///OS deletion
double OsDel(double mu, double beta){
        return fIns(mu,beta,0.0,0.5,-1.0,wghtOne);
}

///EXP insertion
double ExpIns(double mu, double beta){
	return fIns(mu,beta,0.0,1.0,1.0,wghtOne);
}

///EXP deletion
double ExpDel(double mu, double beta){
        return fIns(mu,beta,0.0,1.0,-1.0,wghtOne);
}

///BAR
double BarForBack(double mu, double beta, double parC, double drt){
	return fIns(mu,beta,parC,0.5,drt,wghtBar);
}

//OS
double OsForBack(double mu, double beta, double parC, double drt){
	return fIns(mu,beta,0.0,0.5,drt,wghtOne);
}

///EXP
double ExpForBack(double mu, double beta, double parC, double drt){
	return fIns(mu,beta,0.0,1.0,drt,wghtOne);
}
