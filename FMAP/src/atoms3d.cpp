#include "atoms3d.hpp"

/// init particle positions by assigning them on a cubic grid. Modified from http://www.pages.drexel.edu/~cfa22/msim/codes/mclj_widom.c
void Atoms3D::init_cube(double boxlen){
  int n=1;
  int ix,iy,iz;
  int yoffset=cap,zoffset=cap*2;
  /* Find the lowest perfect cube, n, greater than or equal to the number of particles */
  while ((n*n*n)<cap) n++;

  ix=iy=iz=0;
  for (int i=0;i<cap;i++){
    cords[i] = ((double)ix+0.5)*boxlen/n;
    cords[yoffset+i] = ((double)iy+0.5)*boxlen/n;
    cords[zoffset+i] = ((double)iz+0.5)*boxlen/n;
    ix++;
    if (ix==n){
      ix=0;
      iy++;
      if (iy==n){
	iy=0;
	iz++;
      }
    }
  }
}

///PUTS ATOMS ON A FCC LATTICE,cap <=4*n^3. Modified from SPHRI.F, or ftp://ftp.dl.ac.uk/ccp5/ALLEN_TILDESLEY/F.23
void Atoms3D::init_fcc(double boxlen){
  int n=1;
  int ix,iy,iz,iref;
  int yoffset=cap,zoffset=cap*2;
  while (4*(n*n*n)<cap) n++;
  
  double cell=boxlen/n;
  double cell2=cell/2.0;
  // BUILD THE UNIT cell **

  //  ** SUBLATTICE A **

        cords[0] =  0.0;
        cords[yoffset+0] =  0.0;
        cords[zoffset+0] =  0.0;

   // ** SUBLATTICE B **

        cords[1] =  cell2;
        cords[yoffset+1] =  cell2;
        cords[zoffset+1] =  0.0;

   // ** SUBLATTICE C **

        cords[2] =  0.0;
        cords[yoffset+2] =  cell2;
        cords[zoffset+2] =  cell2;

   // ** SUBLATTICE D **

        cords[3] =  cell2;
        cords[yoffset+3] =  0.0;
        cords[zoffset+3] =  cell2;

  int m=0;
  int cnt=0;
  bool done=false;
  for (ix=0;ix<n && cnt<cap;ix++){
    for (iy=0;iy<n && cnt<cap;iy++){
      for (iz=0;iz<n && cnt<cap;iz++){
	for (iref=0;iref<4 && cnt<cap;iref++){
	  cords[iref+m]=cords[iref]+cell*ix;
	  cords[yoffset+iref+m]=cords[yoffset+iref]+cell*iy;
	  cords[zoffset+iref+m]=cords[zoffset+iref]+cell*iz;
	  cnt++;
	}
	m=m+4;
      }
    }
  }
}

void Atoms3D::init(int n, double xyz[][3]){
  int i;
  if (n>cap){
    cords.resize(n*3);
    cap=n;
  }
  int yoffset=cap,zoffset=cap*2;
  for (i=0;i<cap;i++){
    cords[i]=xyz[i][0];
    cords[yoffset+i]=xyz[i][1];
    cords[zoffset+i]=xyz[i][2]; 
  }
}

void Atoms3D::init(std::vector<Point> xyz){
  int i;
  int n=xyz.size();
  if (n>cap){
    cords.resize(n*3);
    cap=n;
  }
  int yoffset=cap,zoffset=cap*2;
  for (i=0;i<cap;i++){
    cords[i]=xyz[i].x;
    cords[yoffset+i]=xyz[i].y;
    cords[zoffset+i]=xyz[i].z;
  }
}
