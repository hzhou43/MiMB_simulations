#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
using namespace std;

const int M = 5000; 	//max number for array size
int N = 128;  		//number or particles; argument
int boxId[M];		// in which box is a particle? (0 or 1)
const int box = 2;	//number of boxes
int Npbox[2];		//number of particles in each box
double * x = new double[M];	//position arrays
double * y = new double[M];	
double * z = new double[M];
double ucut;
double rcut2;
double beta;
double Utot[box];	//total energy of a box
double Virtot[box];	//total virial of a box
double L[box];		//length of each box
double lmax;
double dVmax;

void initialize();
void ener(int box, double &vir, double &u, double xp, double yp, double zp, int part);
void accept(double factor, int &accept);
void energy(int box, double &vir, double &u);
void particle_disp();
void volume_exch();
void particle_swap();
void read_coor(ostream& traj, int cycle);

int main(int argc, char* argv[]){
	int Ncycle;
	double pdisp, pswap, vexch, temp, vir, u, des, Randomnumber;
	double rho;
	ofstream fout;

	//Input values (arguments)
	Ncycle = 100000;
	N  = atoi(argv[1]);
	rho = atof(argv[3]); //argument
	temp = atof(argv[2]); //argument
	pdisp = N;
	pswap = N;
	vexch = 5;
	lmax  = 0.5;
	dVmax = 0.5;
	L[0] = pow((N/2)/rho, 1./3.);
	L[1] = pow((N/2)/rho, 1./3.);

	beta = 1./temp;
	if (N%2==0){
		Npbox[0] = N/2;
		Npbox[1] = N/2;
	}
	else{
		Npbox[0] = (N+1)/2;
		Npbox[1] = N - Npbox[0] - 1;
	}
	rcut2 = 3.0*3.0;
	ucut  = rcut2*rcut2*rcut2;
	ucut  = 1.0/ucut;
	ucut  = 4.0*ucut*(ucut-1.0);
	
	srand(time(0));
	initialize();
	Utot[0] = Utot[1] = 0;
	Virtot[0] = Virtot[1] = 0;
	for(int i=0; i<2; i++){
		energy(i, vir, u);
		Virtot[i] = vir;
		Utot[i] = u;
	}
	
	ofstream traj;
        traj.open("Trajectory.xyz",ios::out);
	fout.open("output.txt", ios::out);

	for(int cycle=0; cycle<Ncycle; cycle++){
		for(int i=0; i<N; i++){
		des = int(((double)rand()/RAND_MAX)*(pdisp+vexch+pswap));
		if (des < pdisp)
			particle_disp();
		else if (des < pdisp+vexch)
			volume_exch();
		else
			particle_swap();
		}
		fout << cycle << " " <<  Utot[0] << " " << Utot[1];
		fout << " " << pow(L[0],3.0) << " " << pow(L[1], 3.0) << " ";
		fout << Npbox[0] << " " << Npbox[1] << " ";
		fout << Npbox[0]/pow(L[0],3.0) << " " << Npbox[1]/pow(L[1], 3.0);
		fout << "\n";
		
		if(cycle%100==0){
                        read_coor(traj,cycle);
                }
	}
	return 0;
}

void read_coor(ostream& traj, int cycle){
        traj << '\t' << int(N) << "\n\n";
        for(int i=0; i<N; i++){
                if(boxId[i]==0){
			traj << "Ar" << '\t' << x[i]+L[0] << '\t';
                	traj << y[i] << '\t' << z[i] << '\n';
		}
                else{
			traj << "Ne" << '\t' << x[i] << '\t';
                	traj << y[i] << '\t' << z[i] << '\n';
		}
        }
}

void particle_swap(){
	int Laccept=0;
	int add, del, pd;
	double p, uadd, udel, viradd, virdel;
	double xd,yd,zd,xa,ya,za;

	p = (double)rand()/RAND_MAX;

	if(p<0.5){
		add=0;
		del=1;
	}else{
		add=1;
		del=0;
	}

	if(Npbox[del]==1)
		return;

	do{
		pd = int(N*(double)rand()/RAND_MAX);
	}while(boxId[pd]!=del);

	xd = x[pd];
	yd = y[pd];
	zd = z[pd];

	ener(del, virdel, udel, xd, yd, zd, pd);

	xa = L[add]*(double)rand()/RAND_MAX;
	ya = L[add]*(double)rand()/RAND_MAX;
	za = L[add]*(double)rand()/RAND_MAX;

	ener(add, viradd, uadd, xa, ya, za, N);

	accept((double(Npbox[del])*(pow(L[add],3))/((pow(L[del],3))*double(Npbox[add]+1)))*exp(-beta*(uadd-udel)),Laccept);


	if(Laccept){
		Npbox[add] = Npbox[add]+1;
		Npbox[del] = Npbox[del]-1;

		x[pd] = xa;
		y[pd] = ya;
		z[pd] = za;

		boxId[pd] = add;

		
		Utot[add] = Utot[add] + uadd;
		Utot[del] = Utot[del] - udel;

		Virtot[add] = Virtot[add] + viradd;
		Virtot[del] = Virtot[del] - virdel;
	}
}

void volume_exch(){
	int Laccept=0;
	int i;
	double xo[N], yo[N], zo[N], fact[2];
	double xv[N], yv[N], zv[N];
	double un1, un2, virn1, virn2;
	double eno[2], enn[2];
	double vo[2], Lo[2], Ln[2], vn[2];
	double Vn;

	eno[0] = Utot[0];
	eno[1] = Utot[1];
	
	for(i=0; i<2; i++){
		Lo[i] = L[i];
		vo[i] = pow(Lo[i],3.0);;
	}
	Vn = vo[0] + vo[1];
	vn[0] = vo[0] + (2.0*(double)rand()/RAND_MAX-1.0)*dVmax;

	if(vn[0]<0.0)
		return;

	vn[1] = Vn-vn[0];

	if(vn[1]<0.0)
		return;
	
	Ln[0] = pow(vn[0],1./3.);
	Ln[1] = pow(vn[1],1./3.);

	for(int i=0; i<2; i++){
		fact[i] = Ln[i]/Lo[i];
		L[i] = pow(vn[i],1./3.);
	}

	for(int i=0; i<N; i++){
		xo[i] = x[i];
		yo[i] = y[i];
		zo[i] = z[i];

		xv[i] = x[i]*fact[boxId[i]];
		yv[i] = y[i]*fact[boxId[i]];
		zv[i] = z[i]*fact[boxId[i]];
	}

	energy(0, virn1, un1);
	energy(1, virn2, un2);
	
	enn[0] = un1;
	enn[1] = un2;

	accept(exp(-beta*(enn[0]+enn[1]-eno[0]-eno[1]) + double(Npbox[0])*log(vn[0]/vo[0]) + double(Npbox[1])*log(vn[1]/vo[1])), Laccept);


	if(Laccept){

		Utot[0] = un1;
		Virtot[0] = virn1;
	
		Utot[1] = un2;
		Virtot[1] = virn2;
	}
	else{
		L[0] = Lo[0];
		L[1] = Lo[1];
		
		for(int i=0; i<N; i++){
			x[i] = xo[i];
			y[i] = yo[i];
			z[i] = zo[i];
		}
	}
}

void particle_disp(){
	int Laccept=0;	
	int box, p;
	double xd, yd, zd, xo, yo, zo, ran, un, uo;
	double virn, viro;
	
	if(N==0)
		return;
	
	p = int(N*(double)rand()/RAND_MAX);
	box = boxId[p];

	ran = (double)rand()/RAND_MAX;
	xd = x[p] + (2.0*(double)rand()/RAND_MAX-1.0)*lmax;
	yd = y[p] + (2.0*(double)rand()/RAND_MAX-1.0)*lmax;
	zd = z[p] + (2.0*(double)rand()/RAND_MAX-1.0)*lmax;

	if (xd<0.0)
		xd = xd + L[box];
	else if (xd>L[box])
		xd = xd - L[box];
	if (yd<0.0)
                yd = yd + L[box];
        else if (yd>L[box])
                yd = yd - L[box];
	if (zd<0.0)
                zd = zd + L[box];
        else if (zd>L[box])
                zd = zd - L[box];	
	
	xo = x[p];
	yo = y[p];
	zo = z[p];

	ener(box, viro, uo, xo, yo, zo, p);
	ener(box, virn, un, xd, yd, zd, p);

	accept(exp(-beta*(un-uo)), Laccept);

	if(Laccept){
		Utot[box] = Utot[box] + un - uo;
		Virtot[box] = Virtot[box] + virn - viro;
	
		x[p] = xd;
		y[p] = yd;
		z[p] = zd;
	}
}

void energy(int box, double &vir, double &u){
	int i, j;
	double dx, dy, dz, r2, bl, hbl;

	u = 0;
	vir  = 0;

	bl = L[box];
	hbl = bl/2.0;
	
	for (i = 0; i < N-1; i++){
		if (boxId[i]==box){
			for(j=i+1; j<N; j++){
				if (boxId[j]==boxId[i]){
					dx = x[i] - x[j];
					dy = y[i] - y[j];
					dz = z[i] - z[j];
						
					if (dx>hbl)
						dx = dx - bl;
					else if (dx <-hbl)
						dx = dx + bl;
					if (dy>hbl)
                	                        dy = dy - bl;
               		                 else if (dy <-hbl)
                	       	                 dy = dy + bl;
					if (dz>hbl)
                     		                  dz = dz - bl;
                           	     	else if (dz <-hbl)
                                        	dz = dz + bl;
				
					r2 = dx*dx + dy*dy + dz*dz;
			
					if(r2<rcut2){
						r2 = 1.0/r2;
						r2 = r2*r2*r2;
						u = u + 4.0*r2*(r2-1.0)-ucut;
						vir  = vir + 48.0*r2*(r2-0.5);
					}

				}
			}
		}
	}
}

void accept(double factor, int &accept){
	double ran;
	
	ran = (double)rand()/RAND_MAX;
	if (factor > 1.0)
		accept = 1;
	else if ( factor < 0)
		accept = 0;
	else{
		if (ran<factor)
			accept = 1;
		else
			accept = 0;
	}
}

void ener(int box, double &vir, double &u, double xp, double yp, double zp, int part){
	int i;
	double dx, dy, dz, r2, bl, hbl;

	u = 0;
	vir  = 0;
	bl   = L[box];
	hbl  = bl/2.0;

	for ( i=0; i< N; i++){
		if (boxId[i]==box and i!=part){
			dx = x[i] - xp;
			dy = y[i] - yp;
			dz = z[i] - zp;
			
			if (dx > hbl)
				dx = dx - bl;
			else if ( dx < - hbl)
				dx = dx + bl;
			if (dy > hbl)
                                dy = dy - bl;
                        else if ( dy < - hbl)
                                dy = dy + bl;
			if (dz > hbl)
                                dz = dz - bl;
                        else if ( dz < - hbl)
                                dz = dz + bl;

			r2 = dx*dx + dy*dy + dz*dz;

			if ( r2 < rcut2 ){
				r2 = 1.0/r2;
				r2 = r2*r2*r2;
				u = u + 4.0*r2*(r2-1.0)-ucut;
				vir = vir + 48*r2*(r2-0.5);
			}
		}
	}
}

void initialize(){
	int p=0, i=0, j=0, k=0;
		int Nh = N/2;
		int n=0;
		double dx;
		dx = L[0]/pow(Npbox[0],1./3.);

		do{
			n++;
		}while(n*n*n<Nh);

		do{
			while(i<n && p<Nh){
				while(j<n && p<Nh){
					while(k<n && p<Nh){
						x[p]=x[p+Nh-1]=i*dx;
						y[p]=y[p+Nh-1]=j*dx;
						z[p]=z[p+Nh-1]=k*dx;
						boxId[p]=0;
						boxId[p+Nh]=1;
						k++;
						p++;
					}
					k=0;
					j++;
				}
				j=0;
				i++;
			}
			i=0;
		}while(p<Nh);
}
