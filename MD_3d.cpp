#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
#define Npm 1000
#define Pm 500

int copy(double *x_update,double *y_update,double *z_update,double *x,double *y,double *z,int Np){
  int i;
  for(i=0;i<Np;i++){
    x_update[i]=x[i];
    y_update[i]=y[i];
    z_update[i]=z[i];
  }
  return 0;
}


int calc_disp_max(double *disp_max,double *x,double *y,double *z,double *x_update,double *y_update,double *z_update,int Np,int L){
  int i;
  double dx,dy,dz;
  double disp;

  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;

    dy=y[i]-y_update[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;

    dz=z[i]-z_update[i];
    if(dz>0.5*L) dz-=L;
    else if(dz<-0.5*L)dz+=L;
    disp = dx*dx+dy*dy+dz*dz;

    if(disp > *disp_max)
      *disp_max =disp;
  }


  return 0;
}

int com_correction(double *x,double *y,double *z,double *x_corr,double *y_corr,double *z_corr,int Np,double L){
  int i;
  double dx,dy,dz;
  static double x0[Npm],y0[Npm],z0[Npm];
  static bool IsFirst = true;
  if(IsFirst){
    for(i=0;i<Np;i++){
      x0[i]=x[i];
      y0[i]=y[i];
      z0[i]=z[i];
    }
    IsFirst = false;
  }

  for(i=0;i<Np;i++){
    dx=x[i]-x0[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;

    dy=y[i]-y0[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;

    dz=z[i]-z0[i];
    if(dz>0.5*L) dz-=L;
    else if(dz<-0.5*L)dz+=L;

    *x_corr+=dx/Np; //center of mass displacement.x
    *y_corr+=dy/Np;
    *z_corr+=dz/Np;


    x0[i]=x[i];
    y0[i]=y[i];
    z0[i]=z[i];
  }
  return 0;
}

int output(double *x,double *y,double *z,double *r1,int Np,double t,double dt){
  int i;
  static int count=1;
  char filename[128];
  sprintf(filename,"MD_bina2d_%.d_N1000.dat",count);
  ofstream file;
  file.open(filename);

  for(i=0;i<Np;i++)
    file << x[i] << " " << y[i]<< " " << z[i]<<" "<<r1[i] << endl;
  file.close();
  count++;
  return 0;
}

int output_t(double *x,double *y,double *z,double avU, double avK,int Np,double t,double time_stamp,double x_corr,double y_corr,double z_corr)
{
  int i;
  char filename[64];
  char filename2[64];
  FILE *fp;
  FILE *fp2;

  sprintf(filename,"MD_time_coord_N1000.dat");
  fp=fopen(filename,"a+");

  for(i=0;i<Np;i++)
    fprintf(fp,"%f\t%f\t%f\n%f\n",t-time_stamp,x[i]-x_corr,y[i]-y_corr,z[i]-z_corr);
  fclose(fp);

  sprintf(filename2,"MD_time_energy_N1000.dat");
  fp2=fopen(filename2,"a+");
  fprintf(fp2,"%f\t%f\t%f\t%f\t%f\t%f\n%f\n",t-time_stamp, avK,avU ,avK+avU ,x_corr,y_corr,z_corr);
  fclose(fp2);

  return 0;
}

double unif_rand(double left, double right)
{
  return left + (right - left) * rand() / RAND_MAX;
}

double gaussian_rand(void)
{
  static double iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (iset == 0) {
    do {
      v1 = unif_rand(-1, 1);
      v2 = unif_rand(-1, 1);
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 0.50;
    return v2 * fac;
  }
  else {
    iset = 0;
    return gset;
  }
}

int ini_coord_rand(double* x, double* y, double* z,double* a, int Np, double L, double r1, double r2){
  int k, p;
  p = 0;
  for (k = 0; k < Np; k++) {
    x[k] = unif_rand(0, 1) * L;
    y[k] = unif_rand(0, 1) * L;
    z[k] = unif_rand(0, 1) * L;
    if (p == 0) {
      a[k] = r1;
      p = p + 1;
    }
    else {
      a[k] = r2;
      p = p - 1;
    }
  }
  return 0;
}

int ini(double* vx, double* vy,double* vz, int Np) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    vz[j] = 0.0;
  }
  return 0;
}

int f(int i,int M)
{
  int k;

  k=i;

  if(k<0)
    k+=M;
  if(k>=M)
    k-=M;

  return k;
}

void update(double L,int Np,double *x,double *y,double *z,int M,double RCHK,int (*list)[Pm])
{
  int i,j,k;
  int nx,ny,nz;
  int n,l,m;
  double dx,dy,dz,r;

  int (*map)[Npm]=new int[M*M*M][Npm];

  for(k=0;k<M;k++){
    for(j=0;j<M;j++){
      for(i=0;i<M;i++){
     	map[i+M*j+M*M*k][0]=0;
    }
  }
}

  for(i=0;i<Np;i++){
    nx=f((int)(x[i]*M/L),M);
    ny=f((int)(y[i]*M/L),M);
    nz=f((int)(z[i]*M/L),M);//aderess of particle i is determined.

    for(n=nz-1;n<=nz+1;n++){
      for(m=ny-1;m<=ny+1;m++){
        for(l=nx-1;l<=nx+1;l++){
          map[f(l,M)+M*f(m,M)+M*M*f(n,M)][map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0] +1]=i;
          map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0]++;
	//	printf("%d\n", map[f(l,M)+M*f(m,M)][0]);
      }
    }
  }
  }
  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)(x[i]*M/L),M);
    ny = f((int)(y[i]*M/L),M);
    nz = f((int)(z[i]*M/L),M);


    for (k=1; k<=(map[nx+M*ny+M*M*nz][0]); k++){
      j = map[nx+M*ny+M*M*nz][k];
      if(j>i){
        dx =x[i] - x[j];
        dy =y[i] - y[j];
        dz =z[i] - z[j];

        if(dx<-L/2.0)
          dx+=L;
        else if(dx> L/2.0)
          dx-=L;

        if(dy<-L/2.0){
          dy+=L;
        }
        else if(dy> L/2.0)
          dy-=L;
        if(dz<-L/2.0){
          dz+=L;
          }
        else if(dz> L/2.0)
          dz-=L;

        r = dx*dx + dy*dy + dz*dz;

        if(r<RCHK*RCHK){
          list[i][0]++;
          list[i][list[i][0]]=j;
        }
      }
    }
  }
  delete []map;
}

int calc_force_hs(double* x, double* y, double* z,double L, int Np, double* a, double* kx, double* ky,double* kz, double* avU,int (*list)[Pm]) {
  int i, j, k;
  *avU = 0.0;
  double r;
  double t, drU;
  double dx, dy,dz;
  double aij;
  double cut;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kz[k] = 0.0;
  }
  for (i = 0; i < Np; i++){
    for (j = 1; j <=list[i][0]; j++){
      dx = x[list[i][j]] - x[i];
      dy = y[list[i][j]] - y[i];
      dz = z[list[i][j]] - z[i];
      if (dx > (0.5 * L))
        dx -= L;
      if (dx < -(0.5 * L))
        dx += L;
      if (dy > (0.5 * L))
        dy -= L;
      if (dy < -(0.5 * L))
        dy += L;
      if (dz > (0.5 * L))
        dz -= L;
      if (dz < -(0.5 * L))
        dz += L;
      aij = (a[i] + a[list[i][j]]) / 2.0;
      r = sqrt(dx * dx + dy * dy + dz * dz);
      t = r / aij;
      cut = aij;
      if (r < cut) {
        drU = -(1 - t) / aij; //analytical calculation of the 1'st derivative
      }
      else {
        drU = 0.0;
        continue;
      }
      kx[list[i][j]] -= drU * dx / r;
      kx[i] += drU * dx / r;
      ky[list[i][j]] -= drU * dy / r;
      ky[i] += drU * dy / r;
      kz[list[i][j]] -= drU * dz / r;
      kz[i] += drU * dz / r;
      *avU += (1 - t) * (1 - t);
    }
  }
  *avU /= double(Np);
  return 0;
}

int calc_force(double* x, double* y, double* z,double L, int Np, double* a, double* kx, double* ky, double* kz,double* avU,int (*list)[Pm]) {
  int i, j, k;
  *avU = 0.0;
  double r2;
  double w2,w4,w12,drU;
  double dx, dy,dz;
  double aij;
  double cut;
  cut = 3.0;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kz[k] = 0.0;
  }
  for (i = 0; i < Np; i++)
    {
      for (j = 1; j <=list[i][0]; j++)
        {
          dx = x[list[i][j]] - x[i];
          dy = y[list[i][j]] - y[i];
          dz = z[list[i][j]] - z[i];
          if (dx > (0.5 * L))
            dx -= L;
          if (dx < -(0.5 * L))
            dx += L;
          if (dy > (0.5 * L))
            dy -= L;
          if (dy < -(0.5 * L))
            dy += L;
          if (dz > (0.5 * L))
            dz -= L;
          if (dz < -(0.5 * L))
            dz += L;
          aij = (a[list[i][j]] + a[i]) / 2.0;
          r2 = dx * dx + dy * dy + dz * dz; //avoid using sqrt()
          w2 = aij*aij / r2;
          w4=w2*w2;
	  w12=w4*w4*w4;
          if (r2 < cut*cut) {
            drU = (-12.0) * w12 / r2; //analytical calculation of the 1'st derivative / r
          }
          else {
            drU = 0.0;
            continue;
          }
          kx[list[i][j]] -= drU * dx;
          kx[i] += drU * dx;
          ky[list[i][j]] -= drU * dy;
          ky[i] += drU * dy;
          kz[list[i][j]] -= drU * dz;
          kz[i] += drU * dz;

          *avU += w12;
        }
    }
  *avU /= double(Np);
  return 0;
}

int eq_motion_BD(double* x, double* y, double*z,double* vx, double* vy, double* vz,double dt, double* kx, double* ky, double* kz,int Np, double* avK, double Th) {
  double zeta;
  zeta = 1.;
  int k;
  for (k = 0; k < Np; k++) {
    vx[k] += -vx[k] * zeta * dt + kx[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vy[k] += -vy[k] * zeta * dt + ky[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vz[k] += -vz[k] * zeta * dt + kz[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    z[k] += vz[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k];
  }
  *avK = *avK / Np / 2.0;
  return 0;
}

int eq_motion_MD(double* x, double* y,double* z, double* vx, double* vy, double* vz,double dt, double* kx, double* ky, double* kz,int Np, double L, double* a, double* avU, double* avK, int (*list)[Pm]) {
  int k;
  for (k=0; k<Np ; k++){
    vx[k] += 0.5*kx[k] * dt;
    vy[k] += 0.5*ky[k] * dt;
    vz[k] += 0.5*kz[k] * dt;
  } //step1

  calc_force(x,y,z,L,Np,a,kx,ky,kz,&(*avU),list); //step2

  for(k=0; k<Np ; k++){
    vx[k] += 0.5*kx[k] *dt;
    vy[k] += 0.5*ky[k] *dt;
    vz[k] += 0.5*kz[k] *dt;
  } //step3

  for(k=0; k<Np; k++){
    x[k] += vx[k] * dt + 0.5*kx[k] * dt * dt;
    y[k] += vy[k] * dt + 0.5*ky[k] * dt * dt;
    z[k] += vz[k] * dt + 0.5*kz[k] * dt * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k];
  } //step4
  *avK = *avK / Np / 2.0;
  return 0;
}

int p_bound(double* x, double* y, double* z,int Np, double L) {
  int k;
  for (k = 0; k < Np; k++) {
    if (x[k] < 0.0) {
      x[k] = x[k] + L;
    }
    if (x[k] >  L) {
      x[k] = x[k] - L;
    }
    if (y[k] < 0.0) {
      y[k] = y[k] + L;
    }
    if (y[k] >  L) {
      y[k] = y[k] - L;
    }
    if (z[k] < 0.0) {
      z[k] = z[k] + L;
    }
    if (z[k] >  L) {
      z[k] = z[k] - L;
    }
  }
  return 0;
}

int main()
{
  double t, avU=0.0, avK=0.0;
  double x_corr=0.0,y_corr=0.0,z_corr=0.0;
  int i,count=0,count_noupdate=0;
  double sampling_time,time_stamp=0.;
  double disp_max=0.0;
  int Np = 1000;
  double r1=1.0, r2=1.4;
  double dt = 0.002, time_max = 1.e+8; //parameters
  double sampling_time_max=10000.;
  double time_coord=10000.;
  double time_stable = 100. ,time_eq =1000.;
  double Th = 0.90;
  double RCHK=4.5;
  double L = sqrt(double(Np) / 0.8);
  int    M=(int)(L/RCHK);
  //cout << "L=" << L <<" "<< "M="<<M <<endl;

  double* x, * y, * z,* vx, * vy,* vz, * a, * kx, * ky, * kz,*x_update,*y_update, *z_update;
  int (*list)[Pm]=new int[Npm][Pm];

  x = new double[Np];
  y = new double[Np];
  z = new double[Np];
  x_update = new double[Np];
  y_update = new double[Np];
  z_update = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  vz = new double[Np];
  a = new double[Np];
  kx = new double[Np];
  ky = new double[Np];
  kz = new double[Np];

  char filename[128];

  ini_coord_rand(x, y, z,a, Np, L, r1, r2);
  ini(vx, vy, vz,Np);

  //HP
  for (t = 0.; t < time_stable; t += dt) {
    update(L,Np,x,y,z,M,RCHK,list);
    calc_force_hs(x, y,z, L, Np, a, kx, ky,kz, &avU,list);
    eq_motion_BD(x, y, z,vx, vy, vz,dt, kx, ky, kz,Np, &avK, 0.0);
    p_bound(x, y,z, Np, L);
    //    cout << t <<endl;
  }//緩和時間まではHPでしていて、それ以降から平衡時間まではBDであるがポテンシャルをBHHPでやっている。


  //BHHP to equilibrium(BD)
  sampling_time=5.*dt;
  copy(x_update,y_update,z_update,x,y,z,Np);
  for (t = 0; t < time_eq; t += dt) {
    count++;
    update(L,Np,x,y,z,M,RCHK,list);
    calc_force(x, y, z,L, Np, a, kx, ky,kz, &avU,list);
    eq_motion_BD(x, y, z,vx, vy, vz,dt, kx, ky,kz, Np, &avK, Th);
    com_correction(x,y,z,&x_corr,&y_corr,&z_corr,Np, L);
    p_bound(x, y, z,Np, L);
  }


  sprintf(filename,"MD_energy_time_N1000.txt");
  ofstream file;
  file.open(filename);

    //BHHP from equilibrium(MD)
    sampling_time=5.*dt;
    copy(x_update,y_update,z_update,x,y,z,Np);
    for (t = time_eq; t < time_max; t += dt) {
      count++;
      eq_motion_MD(x, y, z, vx, vy, vz, dt, kx, ky, kz, Np, L, a, &avU, &avK, list);
//int eq_motion_MD(double* x, double* y,double* z double* vx, double* vy, double* vz,
//double dt, double* kx, double* ky, double* kz,int Np, double L, double* a, double* avU, double* avK, int (*list)[Pm]) {
      com_correction(x,y,z,&x_corr,&y_corr,&z_corr,Np, L);
      p_bound(x, y, z, Np, L);
      if(count==int(time_coord/dt)){
	output(x,y,z,a,Np,t,dt);
	file <<t<<" "<< avU<<" "<< avK <<endl;
	count=0;
}

    //logarithmic sampling
    if(int(t/dt) == int((sampling_time + time_stamp)/dt)){
      output_t(x,y,z,avU,avK,Np,t,time_stamp,x_corr,y_corr,z_corr);
      sampling_time*=pow(10.,0.1);
      sampling_time=int(sampling_time/dt)*dt;
      if(sampling_time > sampling_time_max/pow(10.,0.1)){
        time_stamp=t;
        sampling_time=5.*dt;
      }
    }

    //auto update
    calc_disp_max(&disp_max,x,y,z,x_update,y_update,z_update,Np,L);
    //  cout <<disp_max <<endl;
    count_noupdate++;
    if(disp_max>0.8*0.8){
      update(L,Np,x,y,z,M,RCHK,list);
      disp_max=0.0;
      copy(x_update,y_update,z_update,x,y,z,Np);
      // cout <<count_noupdate <<endl;
      count_noupdate=0;
    }

  }
  file.close();

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] x_update;
  delete[] y_update;
  delete[] z_update;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] a;
  delete[] kx;
  delete[] ky;
  delete[] kz;
  delete[] list;
  return 0;
}
