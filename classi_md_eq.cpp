# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
using namespace std;
# define Npm 300
# define Pm 500

int a_r1_ex(double* a, int Np, double L, double* r1){
  int k, p;
  p = 0;
  for(k=0;k<Np;k++)
  if (p == 0) {
    a[k] = r1[k];//ここで入レル
    p = p + 1;
  }
  else {
    a[k] = r1[k];
    p = p - 1;
  }
  return 0;
}


//velocityの初期化
int ini(double* vx, double* vy, double* vz,int Np) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    vz[j] = 0.0;
  }
  return 0;
}

//リストブロックの周期境界補正
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

//周期境界条件
int p_bound(double* x, double* y,double* z, int Np, double L) {
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
    if (x[k] >  L) {
      z[k] = z[k] - L;
    }
  }
  return 0;
}

//重心補正
int com_correction(double *x,double *y,double *z, double *x_corr,double *y_corr,double *z_corr,int Np,double L){
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

//ベルリストのアップデート,三次元の領域分割のイメージも後で。
void update(double L,int Np,double* x,double* y,double* z,int M,double RCHK,int (*list)[Pm]){
  int i,j,k;
  int nx,ny,nz;
  int l,m,n;
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
          map[f(l,M)+M*f(m,M)][0]++;
	//	printf("%d\n", map[f(l,M)+M*f(m,M)][0]);
      }
    }
  }
}
}

//Fouce_calclation
int calc_force(double* x, double* y, double*z, double L, int Np, double* a, double* kx, double* ky, double* kz,double* avU,int (*list)[Pm]) {
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


//荷電粒子の運動方程式
int md_motion(double* x, double* y, double* z, double* vx, double* vy, double* vz, double dt, double* kx, double* ky,
  double* kz, int Np, double* a,double L,double* avU,double* avK,int(*list)[Pm]) {

  int k;
  for (k = 0; k < Np; k++) {
    vx[k] += 0.5*kx[k] * dt;
    vy[k] += 0.5*ky[k] * dt;
    vz[k] += 0.5*kz[k] * dt;
  }
  calc_force(x,y,z,L,Np,a,kx,ky,kz,&(*avU),list);
  for(k=0;k<Np;k++){
    vx[k] += 0.5*kx[k] * dt;
    vy[k] += 0.5*ky[k] * dt;
    vz[k] += 0.5*kz[k] * dt;

    x[k] += vx[k] *dt + 0.5*kx[k]*dt*dt;
    y[k] += vy[k] *dt + 0.5*ky[k]*dt*dt;
    z[k] += vz[k] *dt + 0.5*kz[k]*dt*dt;
  }
  *avK = *avK /Np /2.0;
  return 0;
}

//粒子半径や座標を記録
int output(double *x,double *y,double *z,double *r1,int Np,double t,double dt){
  int i;
  static int count=1;
  char filename[128];
  sprintf(filename,"md_bina2d_%.d.dat",count);
  ofstream file;
  file.open(filename);

   for(i=0;i<Np;i++)
    file << x[i] << " " << y[i]<< " " <<z[i]<<" "<< r1[i] << endl;
  file.close();
  count++;
  return 0;
}



int main(){
  double t, avU=0.0, avK=0.0;
  double x_corr=0.0,y_corr=0.0,z_corr=0.0;
  int count=0;
  int i,j,k;
  int Np = 300;
  double time_coord=10000;
  double dt = 0.01, time_max = 1.e+8;
  double RCHK=4.5;
  double L = sqrt(double(Np) / 0.8);
  int    M=(int)(L/RCHK);

  double* x, * y, * z, * vx, * vy, * vz, * kx, * ky, * kz, * a,* r1;

  x = new double[Np];
  y = new double[Np];
  z = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  vz = new double[Np];
  kx = new double[Np];
  ky = new double[Np];
  kz = new double[Np];
  a = new double[Np];
  r1 = new double[Np];


  int (*list)[Pm]=new int[Npm][Pm];
  char filename[64];
  ifstream file;//ファイルの読み取り

  sprintf(filename,"bina2d_9900.dat");
  file.open(filename);

  for(i=0;i<Np;i++){
    file >> x[i] >> y[i] >> z[i] >> r1[i];
  }
  file.close();

  ini(vx, vy, vz, Np);
  a_r1_ex(a,  Np,  L,  r1);

  update(L,Np,x,y,z,M,RCHK,list);//make list

  for(t=dt;t<time_max;t+=dt){
    count++;
    md_motion(x, y,  z,  vx,  vy,  vz,  dt,  kx, ky,
       kz,  Np,  a,  L,  &avU, &avK,list);
       //int md_motion(double* x, double* y, double* z, double* vx, double* vy, double* vz, double dt, double* kx, double* ky,
        // double* kz, int Np, double* a,double L,double* avU,double* avK,int(*list)[Pm])
    com_correction(x,y,z,&x_corr,&y_corr,&z_corr,Np, L);
    p_bound(x, y, z,  Np, L);



    sprintf(filename,"energy_time.txt");
    ofstream file;
    file.open(filename);//ファイルの書き込み

    if(count==int(time_coord/dt)){
      output(x,y,z,a,Np,t,dt);
      file <<t<<" "<< avU<<" "<< avK <<endl;
      count=0;
    }
  }





  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] kx;
  delete[] ky;
  delete[] kz;
  delete[] r1;
  delete[] list;
  delete[] a;
  return 0;
}
