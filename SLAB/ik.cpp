#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <string>
#include <omp.h>

using namespace std;

int findBin(double x, int M){
	return (int)(floor(x));	
}

int readCoor(ifstream &fin, int N, double coor[][4]){
	int i=0;
	while(i<N){
		fin >> coor[i][0] >> coor[i][1] >> coor[i][2] >> coor[i][3];
		i++;
	}
	return i;
}

void cenCoor(int N, double coor[][4]){
	double xcm, ycm, zcm;
	double tol=0.1;
	xcm=ycm=zcm=0;
	//Calculate CM and translate at the center of the box
	for(int p=0; p<N; p++){
		zcm += coor[p][2];
	}
	zcm = zcm/N;
	while(abs(zcm)>tol){
		for (int p=0; p<N; p++){
			coor[p][2] -= zcm;
			if (coor[p][2] < -75){
				coor[p][2] += 150;
			}
			if (coor[p][2] > 75){
				coor[p][2] -= 150;
			}
		}
		zcm=0;
		for (int p=0; p<N; p++){
			zcm += coor[p][2];
		}
		zcm = zcm/N;
	}
}

void denHist(int N, double coor[][4], double Lz, double dz, int M, double histT[], double histA[], double histB[]){
	for(int p=0; p<N; p++){
		int b=findBin((coor[p][2]+Lz/2.0)/dz,M);
		histT[b] += 1;
		if(coor[p][3]==0){
			histA[b] += 1;
		}else{ 
			histB[b] += 1;
		}
	}
}

double calcU(double eps1, double eps2, double eps12, double rij, double rij2, double pair){
	int type=(int)pair;
	double u;
	double r6 = rij2 * rij2 * rij2;
	double r12 = r6 * r6;
	switch(type){
		case 1:
			u=24 * eps12 * (1.0/rij) * (2*(1.0/r12) - (1.0/r6));	
			break;
		case 0:
			u=24 * eps1 * (1.0/rij) * (2*(1.0/r12) - (1.0/r6));
			break;
		case 2:
			u=24 * eps2 * (1.0/rij) * (2*(1.0/r12) - (1.0/r6));
			break;
		default:
			cerr <<"Pair error. Expected: 0,1,2"<<" Got: "<<type<<"\n";
			exit(-1);	
	}
	return u;
}

double pbcR(double r, double l){
	if (r>l/2.0){
		r = r - l;
	}else if(r<-l/2.0){
		r = r + l;
	}
	return r;
}

void press(int M, double pn[], double pt[], int bi, int bj, double zi, double zj, double xr, double yr,double zr, double rij, double devU, double Lz){
	int nij;
	if( bi <= bj and abs(zi-zj)<Lz/2.0 ){
		nij = (abs(bj-bi)-1) + 2;
		for(int pb=bi; pb<=bj; pb++){
			double mypn=((zr*zr)/rij)*devU*(1.0/nij);
			double mypb=(0.5*(xr*xr+yr*yr)/rij)*devU*(1.0/nij);
			pn[pb] += mypn;
			pt[pb] += mypb;
		}
	}
	else if( bi>bj and abs(zi-zj)<Lz/2.0){
		nij = (abs(bi-bj) -1) + 2;
		for(int pb=bj; pb<=bi; pb++){
			double mypn=((zr*zr)/rij)*devU*(1.0/nij);
			double mypb=(0.5*(xr*xr+yr*yr)/rij)*devU*(1.0/nij);
			pn[pb] += mypn;
			pt[pb] += mypb;
		}
	}//Account MIC
	else if( bi>bj and abs(zi-zj)>Lz/2.0 ){
		nij = abs(M-bi+bj);
		for(int pbi=bi;pbi<M;pbi++){
			double mypn=((zr*zr)/rij)*devU*(1.0/nij);
			double mypb=(0.5*(xr*xr+yr*yr)/rij)*devU*(1.0/nij);
			pn[pbi] += mypn;
			pt[pbi] += mypb;
		}
		for(int pbj=0; pbj<=bj; pbj++){
			double mypn=((zr*zr)/rij)*devU*(1.0/nij);
                        double mypb=(0.5*(xr*xr+yr*yr)/rij)*devU*(1.0/nij);
			pn[pbj] += mypn;
			pt[pbj] += mypb;
		}
	}else if( (bi<=bj) and abs(zj-zi)>Lz/2.0 ){
		nij = abs(M-bj+bi);
		for(int pbi=0; pbi<=bi; pbi++){
			double mypn=((zr*zr)/rij)*devU*(1.0/nij);
                        double mypb=(0.5*(xr*xr+yr*yr)/rij)*devU*(1.0/nij);
			pn[pbi] += mypn;
			pt[pbi] += mypb;
		}
		for(int pbj=bj; pbj<M; pbj++){
			double mypn=((zr*zr)/rij)*devU*(1.0/nij);
                        double mypb=(0.5*(xr*xr+yr*yr)/rij)*devU*(1.0/nij);
			pn[pbj] += mypn;
			pt[pbj] += mypb;
		}
	}
	else{
		cerr <<"Never found\n";
		exit(-1);
	}
}

int main(int argc, char* argv[]){
	int N=20000; //number of particles
	double coor[N][4]; //system coordinates matrix, reading from Tcoor.txt
	double eps12=0.30, T=0.65;
	double Lx, Ly, Lz, S;
	int num=9,traj=1,trj=1, M=1;
	double rcut, rcut2, dz;

	eps12=atof(argv[1]);
	M=atoi(argv[2]);
	trj=atoi(argv[3]);
	double histT[M], histA[M], histB[M], pn[M], pt[M]; //histogram arrays for densities and pressures (pn-normal, pt-transverse)
	
	//initialize arrays
	for(int i=0; i<M; i++){
		histT[i]=0;
		histA[i]=0;
		histB[i]=0;
		pt[i]=0;
		pn[i]=0;
	}

	ifstream fin;
	ofstream fout;
	fin.open("Tcoor.txt"); //coordinates
	if(!fin){
		cerr << "Attempt to open file failed!\n";
		exit(1);
	}
	
	Lz=150;    //box length longitudinal
	Lx=Ly=30;  //box cross section
	dz = Lz/M; // bin size along Ly
	rcut=3.0; 
	rcut2 = rcut*rcut;

	for(int i=0; i<N; i++){
		coor[i][0]=0;
		coor[i][1]=0;
		coor[i][2]=0;
		coor[i][3]=0;
	}
	num=0;	
	while(!fin.eof()){
		fin >> traj;
		int nr=readCoor(fin, N, coor);
		if (nr!=N){
			cerr <<"Corrupted coor file\n";
			cerr <<"Expected: "<<N<<" Got: "<<nr<<"\n";
			exit(-1);
		}
		if( traj>trj ){
			cenCoor(N,coor);
			//Density Histrogram
			denHist(N,coor,Lz,dz,M,histT,histA,histB);
			num = num + 1;
			//Calculate pressure
			for(int i=0; i<N-1; i++){
				double xi, yi, zi;
				int bi;
				xi = coor[i][0];
				yi = coor[i][1];
				zi = coor[i][2];
				bi=findBin((zi+Lz/2.0)/dz,M);
				for(int j=i+1; j<N; j++){
					double xj, yj, zj;
					double xr, yr, zr;
					int bj;
					xj = coor[j][0];
					yj = coor[j][1];
					zj = coor[j][2];
					xr = xi-xj;
					yr = yi-yj;
					zr = zi-zj;
				
					//MIC
					xr=pbcR(xr,Lx);
					yr=pbcR(yr,Ly);
					zr=pbcR(zr,Lz);
					double rij2 = xr*xr + yr*yr + zr*zr;
					if(rij2 <= rcut2){
						double rij = sqrt(rij2);
						double devU=calcU(1.0, 0.9, eps12, rij, rij2, coor[i][3]+coor[j][3]);

						//Find the bins that yi and yj belond to
						bj=findBin((zj+Lz/2.0)/dz,M);
						press(M,pn,pt,bi,bj,zi,zj,xr,yr,zr,rij,devU,Lz);
					}

				}

			}
		}
		cout << traj << endl;
		traj++;
	}
	S = Lx*Ly;                   //cross section of simulation box
	fout.open("histograms.txt"); //Output file, last column is the quantity of interest (pn-pt -> surface tension)
	for(int k=0; k<M; k++){
			fout << k*dz-Lz/2.0 << ' ' << (histT[k]/num)/(dz*S) << ' ' <<  (histA[k]/num)/(dz*S) << ' ' << (histB[k]/num)/(dz*S) << ' ' << (1/(S*dz))*(1.0*pn[k])*(1.0/num)-(1/(S*dz))*(1.0*pt[k])*(1.0/num) << '\n';
	}

	return 0;
}
