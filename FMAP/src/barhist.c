/**
 *    @file    barhist.c
 *    @author  Sanbo Qin
 *    @brief   Histgram for BAR.
 */
#include<stdio.h>
#include<math.h>
#include<stddef.h>
#include<stdlib.h>
#include<malloc.h>
#include<assert.h>
#include<limits.h>
#include<stdbool.h>

#include "bar.h"
#include "barhist.h"

struct BarHistST {
	///direction: forward 1.0, backward -1.0, i.e value save as mu or -mu
	double drt;
	///number of bins
	int n;
	///low boundary
	double low;
        ///bin width
	double bw;
	///container
	double *h;
};

BarHist BarHistInit(double drt, int n, double low, double bw){
	if (n<=0 ||n >INT_MAX){
		fprintf(stderr,"ERROR at %s:%d, nbin: %d < %d\n",__FILE__,__LINE__,n,0);
		return NULL;
	}
	if (bw<=0.0){
		fprintf(stderr,"ERROR at %s:%d, binwith: %f < %f\n",__FILE__,__LINE__,bw,0.0);
		return NULL;
	}
	BarHist b;
	if(!(b = malloc(sizeof(struct BarHistST)))) {
                return NULL;
        }
	if(!(b->h = malloc(sizeof(double)*n))) {
                return NULL;
        }
	b->drt=drt;
	b->n=n;
	b->low=low;
	b->bw=bw;
	return b;
}

void BarHistFree(BarHist *b){
	assert(b && *b);
	free((*b)->h);
	free(*b);
	(*b)=NULL;
}

void BarHistZero(BarHist b){
	int i;
	if (b){
		for (i=0;i<b->n;i++){
			b->h[i]=0.0;
		}
	}
}

void BarHistUpdate(BarHist b, double value){
	BarHistUpdateWght(b,value,1.0);
}

void BarHistUpdateWght(BarHist b, double value, double w){
	double val=b->drt*value;
	if (val<b->low){
		fprintf(stderr,"ERROR at %s:%d, val: %f < lowlim %f\n",__FILE__,__LINE__,val,b->low);
		exit(EXIT_FAILURE);
	}
	assert(val>=b->low);
	double di=round((val-b->low)/b->bw);
	int i=(di<(double)(b->n))?(int)(di):(b->n-1);
	if (!(i>=0)){
		fprintf(stderr,"ERROR at %s:%d, i: %d < %d\n",__FILE__,__LINE__,i,0);
		fprintf(stderr,"ERROR at %s:%d, lowlim: %f bw: %f\n",__FILE__,__LINE__,b->low,b->bw);
	}
	assert(i>=0);
	assert(i<b->n);
	b->h[i]+=w;
}

double BarHistSum(BarHist b, double beta, double parC, double (*barfwbw)(double, double, double, double), double *nsum){
	double bsum=0.0;
	*nsum=0.0;
	int i;
	if (b){
		for (i=0;i<b->n;i++){
			//for valid data, b->h[i] should >=1.0;
			if ((b->h[i])>0.5){
				double mu=(b->low+i*b->bw)*b->drt;
				double cnt=b->h[i];
				bsum+=cnt*barfwbw(mu,beta,parC,b->drt);
				(*nsum)+=cnt;
			}
		}
	}else{
		bsum=0.0;
	}
	return bsum;
}

void BarHistSave(BarHist b, char* fname, bool txt){
	FILE* fp;
	if (txt){
		fp=fopen(fname,"w");
		int i;
		fprintf(fp,"#%19lf%20d%20lf%20lf\n",b->drt,b->n,b->low,b->bw);
		for (i=0;i<b->n;i++){
			double mu=(b->low+i*b->bw)*b->drt;
			double cnt=b->h[i];
			fprintf(fp,"%20lE%20lE\n",mu,cnt);
		}
	}else{
		fp=fopen(fname,"wb");
		double bhead[4];
		bhead[0]=b->drt;bhead[1]=b->n;bhead[2]=b->low;bhead[3]=b->bw;
		fwrite(bhead,sizeof(double),4,fp);
		fwrite(b->h,sizeof(double),b->n,fp);
	}
	if (fp) fclose(fp);
}

#define MAXLENLINE 255
BarHist BarHistLoad(char* fname, bool txt){
	FILE* fp;
	BarHist b=NULL;
	if (txt){
		fp=fopen(fname,"r");
		char buf[MAXLENLINE];
		int n;
		double drt,low,bw,mu,val;
        	if(fgets(buf,MAXLENLINE,fp)!= NULL){
			sscanf(buf+5,"%lf%d%lf%lf\n",&drt,&n,&low,&bw);
			b=BarHistInit(drt,n,low,bw);
		}
		int cnt=0;
		while (fgets(buf,MAXLENLINE,fp)!= NULL){
			sscanf(buf,"%lf%lf\n",&mu,&val);
			b->h[cnt]=val;
			cnt++;
		}
		assert(n==cnt);
	}else{
		fp=fopen(fname,"rb");
		double bhead[4];
		fread(bhead,sizeof(double),4,fp);
		int n=(int)(bhead[1]);
		assert(n>=0);
		b=BarHistInit(bhead[0],n,bhead[2],bhead[3]);
		int ni=fread(b->h,sizeof(double),n,fp);
		assert(ni==n);
	}
	return b;
}

double osBarHist2(BarHist forw, BarHist back, double beta){
	double w0,n0,w1,n1;
	w0=BarHistSum(forw,beta,0.0,OsForBack,&n0);
	w1=BarHistSum(back,beta,0.0,OsForBack,&n1);
	return log(w1/n1)-log(w0/n0);
}

double expBarHist2(BarHist forw, double beta){
	double w0,n0;
	w0=BarHistSum(forw,beta,0.0,ExpForBack,&n0);
	return -log(w0/n0);
}

double barBarHist2(BarHist forw, BarHist back, double beta, double parC){
	double w0,n0,w1,n1;
	w0=BarHistSum(forw,beta,parC,BarForBack,&n0);
        w1=BarHistSum(back,beta,parC,BarForBack,&n1);
	return log(w1/n1)-log(w0/n0);
}

double barIter(BarHist forw, BarHist back, double beta, double parC){
	double bar,nextC;
	double w0,n0,w1,n1;
	w0=BarHistSum(forw,beta,parC,BarForBack,&n0);
        w1=BarHistSum(back,beta,parC,BarForBack,&n1);
	bar=log(w1/n1)-log(w0/n0);
	nextC=(bar+log(n1/n0))/beta;
	return nextC;
}

void BarHistOpt(BarHist forw, BarHist back, double beta, double tol,double boe[6]){
	BarHistOptMore(forw,back,beta,tol,6,boe);
}

void BarHistOptMore(BarHist forw, BarHist back, double beta, double tol, int n, double boe[]){
	double os_f,os_b;
	double ex,os,bar,hf,hb;
	double n0,n1;
	assert(forw->drt>0.0 && back->drt<0.0);
	
	ex=expBarHist2(forw,beta);

        os_f=BarHistSum(forw,beta,0.0,OsForBack,&n0);
	os_b=BarHistSum(back,beta,0.0,OsForBack,&n1);
	os=log(os_b/n1)-log(os_f/n0);
	hf=-log(os_f/n0);
	hb=-log(os_b/n1);
	
	//double parC=(os+log(n1/n0))/beta;
	double parC=0.0;
	double error=1.0;
	int iter=0;
	double converged=1.0;
	while (error>tol){
		double nextC=barIter(forw,back,beta,parC);
		error=fabs(nextC-parC);
		parC=nextC;
		iter++;
		if (iter>20){
			fprintf(stderr,"Coundn't converger after %d iterations, error: %E > tol: %E\n",iter,error,tol);
			converged=-1.0;
			break;
		}
	}
	bar=barBarHist2(forw,back,beta,parC);
	//printf("%16d%16E%16E%16E%16d\n",(int)(n0),bar/beta,os/beta,ex/beta,(int)(n1));
	boe[0]=bar;
	boe[1]=os;
	boe[2]=ex;
	boe[3]=n0;
	boe[4]=n1;
	boe[5]=converged;
	if (n==8){
		boe[6]=hf;
		boe[7]=hb;
	}
}

void BarHistOptLoad(char* fnforw, char* fnback, bool txt, double beta, double tol, double boe[6]){
	BarHist forw=BarHistLoad(fnforw,txt);
	BarHist back=BarHistLoad(fnback,txt);
	BarHistOpt(forw,back,beta,tol,boe);
	BarHistFree(&forw);
	BarHistFree(&back);
}

void BoeRep(FILE *fp, int n, double res[],double beta){
        double bar,os,ex,n0,n1,converged;
        bar=res[0];os=res[1];ex=res[2];n0=res[3];n1=res[4];converged=res[5];
	if ((n0<1e15) && (n1<1e15)){
        	fprintf(fp,"%16.0f%16E%16E%16E%16.0f%16.0f",n0,bar/beta,os/beta,ex/beta,n1,converged);
	}else{
		fprintf(fp,"%16e%16E%16E%16E%16e%16.0f",n0,bar/beta,os/beta,ex/beta,n1,converged);
	}
        if (n==8){
        	double hf,hb;
		hf=res[6];hb=res[7];
        	fprintf(fp,"%16E%16E",hf/beta,hb/beta);
        }
        fprintf(fp,"\n");
}
