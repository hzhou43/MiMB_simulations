#include <math.h>
#include "grd.h"

/*One dimesion mapping into the PBC box*/
int PBC(const int x, const int mx){
        int px;
        if (x>=0&&x<mx){
                return x;
        }else{
                px=x%mx;
                return px>=0?px:px+mx;
        }
}

void floorfrac(int n, double x[], int xi[], double xf[]){
  int i;
  double fxi;
  for (i=0;i<n;i++){
    xi[i]=(int)floor(x[i]);
    xf[i]=x[i]-xi[i];
  }
}

void dst1(const int n, const double f[], double d[],const int shift, const double dx){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		d[i]=(shift-f[i])*dx;
	}
}

void dst2init(const int n, const double xd[],double dst2[]){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		dst2[i]=xd[i]*xd[i];
	}
}

void dst2add(const int n, const double xd[],double dst2[]){
        #pragma omp for schedule(static) nowait
        for (int i=0;i<n;i++){
                dst2[i]+=xd[i]*xd[i];
	}
}

void dst3(const int n, const double xf[], const double yf[], const double zf[], double dst2[], double tmp[],\
		const int xshift, const int yshift, const int zshift, const double dx){
	dst1(n,xf,tmp,xshift,dx);
	dst2init(n,tmp,dst2);
	dst1(n,yf,tmp,yshift,dx);
	dst2add(n,tmp,dst2);
	dst1(n,zf,tmp,zshift,dx);
	dst2add(n,tmp,dst2);
}

void pbc1(const int n, const int x[], int px[], const int lc, const int shift){
	#pragma omp for schedule(static) nowait	
	for (int i=0;i<n;i++){
		px[i]=PBC(x[i]+shift,lc);
#ifdef DEBUG
        //fprintf(stderr,"%s\t:%d\t:lc:%d\tx[i]%dpx[i]%d\n",__FILE__,__LINE__,lc,x[i],px[i]);
	assert(px[i]<lc);
	assert(px[i]>=0);
#endif
	}
}

void dstflt(const int n, const double d2[], double d2f[], const double low2){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		d2f[i]=(d2[i]>low2)?d2[i]:low2;
	}
}

void dst2val(const int n, const double d2[], const double rcut2,const double coef1[], const double coef2[],\
	     double v[], double (*rdf)(double,double)){
	#pragma omp for schedule(static) nowait
	for (int i=0;i<n;i++){
		v[i]=(d2[i]>rcut2)?0.0:(coef1[i]*rdf(d2[i],coef2[i]));
	}
}

void val2grd(const int n, const double v[], const int px[], const int py[], const int pz[],\
	     const int l, fftw_real* grid){
	#pragma omp for schedule(static) nowait
	int N=2*(l/2+1);
	for (int i=0;i<n;i++){
		int ijk=(px[i]*l+py[i])*N+pz[i];
		#pragma omp atomic
		grid[ijk]+=v[i];
	}
}

double estmin(const int i, const double dx){
	if (i<0) return i*dx;
	if (i>1) return (i-1)*dx;
	return 0.0; 
}

double findmaxval(const int n, const double val[]){
	int i;
	if (n==1){ 
		return val[0];
	}else{
		double mx=val[0];
		//#pragma omp parallel for schedule(static) reduction(max : mx)
		for (i=0;i<n;i++){
			if (val[i]>mx){
				mx=val[i];
			}
		}
		return mx;
	}
}

//xi=floor(x), xf=x-xi, xf in [0,1) fill between [ceil(rcut/dx)...ceil(rcut/dx)+1]
void grdR(const int l, fftw_real* grid, \
         const int n, int xi[], int yi[], int zi[], double xf[], double yf[], double zf[], \
	 double coef1[],double coef2[],\
         const double rup, const double rlow, const double dx, double (*rdf)(double,double)){
	int i,j,k;
	double d2[n],d2f[n];
	int px[n],py[n],pz[n];
	int start, end;
	double rcut2=rup*rup;
	double low2=rlow*rlow;
	double estx,esty,estz; //for early skip
	end=(int)(ceil(rup/dx))+1;
	start=-end+1;
#ifdef DEBUG
        printf("%s\t:%d\t:Start:%d\tEnd:%d\n",__FILE__,__LINE__,start,end);
	printf("rup:%f rlow:%f dx:%f\n",rup,rlow,dx);
	printf("%i %i %i %f %f %f\n",xi[n-1],yi[n-1],zi[n-1],xf[n-1],yf[n-1],zf[n-1]);
#endif
	for (i=start;i<=end;i++){
		#pragma omp parallel
		pbc1(n,xi,px,l,i);
		estx=estmin(i,dx);
		for(j=start;j<=end;j++){
			#pragma omp parallel
			pbc1(n,yi,py,l,j);
			esty=estmin(j,dx);
			if ((estx*estx+esty*esty)>rcut2) continue;
			for (k=start;k<=end;k++){
				estz=estmin(k,dx);
				if ((estx*estx+esty*esty+estz*estz)>rcut2) continue;
				#pragma omp parallel
				{
				pbc1(n,zi,pz,l,k);
				dst3(n,xf,yf,zf,d2,d2f,i,j,k,dx);
				dstflt(n,d2,d2f,low2);
				dst2val(n,d2f,rcut2,coef1,coef2,d2,rdf);
				val2grd(n,d2,px,py,pz,l,grid);
				}
			}
		}
	}
}
