# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
using namespace std;
# define Npm 300
int output_M(double *Lx,double *Ly,double *Lz,int Np,double t,int ex_num,double L){
  int i;
  char filename[64];
  FILE *fp;

  sprintf(filename,"orbital_magnetic_moment_Np_%d_ex_num_%d.dat",Np,ex_num);
  fp=fopen(filename,"a+");
  *Lx=*Lx/(L*L*L);
  *Ly=*Ly/(L*L*L);
  *Lz=*Lz/(L*L*L);
  fprintf(fp,"%f\t%f\t%f\t%f\n",t,*Lx,*Ly,*Lz);
    fclose(fp);
 return 0;

}
int output(double *x,double *y,double *z,double *vx,double *vy,double *vz,double t,int ex_num,int Np){
  int i;
  char filename[64];
  FILE *fp;
 char filename2[64];
  FILE *fp2;
char filename3[64];
  FILE *fp3;
char filename4[64];
  FILE *fp4;
  double ave_Vx=0.,ave_Vy=0.,ave_Vz=0.;
  double sigma_Vx=0.,sigma_Vy=0.,sigma_Vz=0.;

 sprintf(filename,"trajectory_one_particle_ex_num_%d.dat",ex_num);
  fp=fopen(filename,"a+");
 sprintf(filename2,"velocity_one_particle_ex_num_%d.dat",ex_num);
  fp2=fopen(filename2,"a+");
sprintf(filename3,"velocity_average_Np_300_ex_num_%d.dat",ex_num);
  fp3=fopen(filename3,"a+");
sprintf(filename4,"velocity_variance_Np_300_ex_num_%d.dat",ex_num);
  fp4=fopen(filename4,"a+");

  fprintf(fp,"%f\t%f\t%f\t%f\n",t,x[0],y[0],z[0]);
 fprintf(fp2,"%f\t%f\t%f\t%f\n",t,vx[0],vy[0],vz[0]);
 for(i=0;i<Np;i++){
   ave_Vx+=vx[i]/Np;
   ave_Vy+=vy[i]/Np;
   ave_Vz+=vz[i]/Np;
}
 fprintf(fp3,"%f\t%f\t%f\t%f\n",t,ave_Vx,ave_Vy,ave_Vz);
 for(i=0;i<Np;i++){
   sigma_Vx+=(vx[i]-ave_Vx)*(vx[i]-ave_Vx)/Np;
   sigma_Vy+=(vy[i]-ave_Vy)*(vy[i]-ave_Vy)/Np;
   sigma_Vz+=(vz[i]-ave_Vz)*(vz[i]-ave_Vz)/Np;
}

 sigma_Vx=pow(sigma_Vx,0.5);
 sigma_Vy=pow(sigma_Vy,0.5);
 sigma_Vz=pow(sigma_Vz,0.5);
 fprintf(fp4,"%f\t%f\t%f\t%f\n",t,sigma_Vx,sigma_Vy,sigma_Vz);
  fclose(fp);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
 return 0;
}
int all_particle_coord_see(double *x,double *y,double *z,int Np,int ex_num,double* avK,double L){
  int i;
  char filename[64];
  FILE *fp;

  sprintf(filename,"all_particle_coord_initial_state_ex_num_%d.dat",ex_num);
  fp=fopen(filename,"a+");
  double ave_x=0,ave_y=0,ave_z=0;
  for(i=0;i<Np;i++){
    ave_x+=x[i]/Np;
    ave_y+=y[i]/Np;
    ave_z+=z[i]/Np;
}
  for(i=0;i<Np;i++)
    fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",x[i],y[i],z[i],*avK,L,ave_x,ave_y,ave_z);
  fclose(fp);

  return 0;

}

int output_t(double *x,double *y,double *z,int Np,double t,double time_stamp,int ex_num)
{
  int i;
  char filename[64];
  FILE *fp;

  sprintf(filename,"time_coord_ex_num_%d.dat",ex_num);
  fp=fopen(filename,"a+");

  for(i=0;i<Np;i++)
    fprintf(fp,"%f\t%f\t%f\t%f\n",t-time_stamp,x[i],y[i],z[i]);
  fclose(fp);

  return 0;
}

double seed_change(int ex_num){
  srand((unsigned) time(NULL)+ex_num);
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
return 0;
}

  int calc_force_hs(double* x, double* y,double* z, double L, int Np, double* a, double* kx, double* ky, double* kz,double* avU){
    int i, j, k;
    *avU = 0.0;
    double r;
    double t, f;
    double dx, dy,dz;
    double aij;
    double cut;
    for (k = 0; k < Np; k++) {
      kx[k] = 0.0;
      ky[k] = 0.0;
      kz[k] = 0.0;
    }
    for (i = 0; i < Np; i++){
       for (j = i + 1; j < Np; j++){
         dx = x[i] - x[j];
         dy = y[i] - y[j];
         dz = z[i] - z[j];
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
          aij = (a[i] + a[j]) / 2.0;
          r = sqrt(dx * dx + dy * dy + dz * dz);
          t = r / aij;
          cut = aij;
          if (r < cut) {
  	  f = -(1 - t) / aij; //analytical calculation of the 1'st derivative
          }
          else {
  	  f = 0.0;
  	  continue;
          }
          kx[i] -= f * dx / r;
          kx[j] += f * dx / r;
          ky[i] -= f * dy / r;
          ky[j] += f * dy / r;
          kz[i] -= f * dz / r;
          kz[j] += f * dz / r;
          *avU += (1 - t) * (1 - t);
        }
    }
    *avU /= double(Np);
    return 0;
  }

int ini_box_rand(double* x, double* y, double* z,double* a, int Np, double L, double r1, double r2){
  int k, p;
  p = 0;
  for (k = 0; k < Np; k++) {
    x[k] = unif_rand(-0.5, 0.5) * L;
    y[k] = unif_rand(-0.5, 0.5) * L;
    z[k] = unif_rand(-0.5, 0.5) * L;
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

int ini(double* vx, double* vy, double* vz,int Np) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    vz[j] = 0.0;
  }
  return 0;
}

int ini_force(double* kx,double* ky,double* kz,int Np){
int j;
for(j=0;j<Np;j++){
  kx[j]=0.0;
  ky[j]=0.0;
  kz[j]=0.0;
}
  return 0;
}


int eq_motion(double* x, double* y,double* z, double* vx, double* vy,double* vz, double* kx,double* ky,double* kz,double dt, int Np, double* avK,double* Lx,double* Ly,double* Lz, double Th) {
  double Z;
  double B;
  double μ;
  Z =1.;
  B = 0;
  μ = 1.;//無次元化した際に出てくる係数μ*
  int k;
  double Vx;//vxだけ新しいものを用いられるのを防ぐ

  *Lx=0;
  *Ly=0;
  *Lz=0;
  for (k = 0; k < Np; k++) {
    Vx=vx[k];
    vx[k] += -vx[k] * Z * dt + kx[k]*dt + vy[k] * B * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();//磁場による項を追加
    vy[k] += -vy[k] * Z * dt + ky[k]*dt - Vx * B * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();
    vz[k] += -vz[k] * Z * dt + kz[k]*dt + sqrt(2. * Z * Th * dt) * gaussian_rand();
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    z[k] += vz[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k];
    *Lx +=  μ * (y[k] * vz[k] - z[k] * vy[k]);
    *Ly +=  μ * (z[k] * vx[k] - x[k] * vz[k]);
    *Lz +=  μ * (x[k] * vy[k] - y[k] * vx[k]);
  }
  *avK = *avK / Np / 2.0;
  return 0;
}
int p_bound(double* x, double* y,double* z, int Np, double L) {
  int k;
  for (k = 0; k < Np; k++) {
    if (x[k] < -0.5*L) {
      x[k] = x[k] + L;
    }
    if (x[k] > 0.5*L) {
      x[k] = x[k] - L;
    }
    if (y[k] < -0.5*L) {
      y[k] = y[k] + L;
    }
    if (y[k] > 0.5*L) {
      y[k] = y[k] - L;
    }
      if (z[k] < -0.5*L) {
        z[k] = z[k] + L;
      }
      if (z[k] > 0.5*L) {
        z[k] = z[k] - L;
    }
  }
  return 0;
}

int  non_equili_eq_motion(double* x, double* y,double* z, double* vx, double* vy,double* vz, double* kx,double* ky,double* kz,double dt, int Np, double* avK,double* Lx,double* Ly,double* Lz,double L) {
  double Z;
  double B;
  double μ;
  Z =1.;
  B = 3;
  μ = 1.;//無次元化した際に出てくる係数μ*
  int k;
  double Vx;//vxだけ新しいものを用いられるのを防ぐ
  double Th_x,Th_y,Th_z;//平衡崩し
  *Lx=0;
  *Ly=0;
  *Lz=0;
  for (k = 0; k < Np; k++) {
    Vx=vx[k];
    if(x[k]<0){
      Th_x = (8/L) * x[k] + 5;}
    if(y[k]<0){
      Th_y = (8/L) * y[k] + 5;}
    if(z[k]<0){
      Th_z = (8/L) * z[k] + 5;}
    if(x[k]>=0){
      Th_x = -(8/L) * x[k] + 5 ;}
    if(y[k]>=0){
      Th_y = -(8/L) * y[k] + 5;}
    if(z[k]>=0){
      Th_z = -(8/L) * z[k] + 5;}
    vx[k] += -vx[k] * Z * dt + kx[k]*dt + vy[k] * B * dt + sqrt(2. * Z * Th_x * dt) * gaussian_rand();//磁場による項を追加
    vy[k] += -vy[k] * Z * dt + ky[k]*dt - Vx * B * dt + sqrt(2. * Z * Th_y * dt) * gaussian_rand();
    vz[k] += -vz[k] * Z * dt + kz[k]*dt + sqrt(2. * Z * Th_z * dt) * gaussian_rand();
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    z[k] += vz[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k];
    *Lx +=  μ * (y[k] * vz[k] - z[k] * vy[k]);
    *Ly +=  μ * (z[k] * vx[k] - x[k] * vz[k]);
    *Lz +=  μ * (x[k] * vy[k] - y[k] * vx[k]);
  }
  *avK = *avK / Np / 2.0;
  return 0;
}

int coord_average(double*x,double*y,double*z,int Np,double t,int ex_num){
  int i;
  double av_x=0,av_y=0,av_z=0;
  char filename[64];
  FILE *fp;
  for(i=0;i<Np;i++){
    av_x+=x[i]/Np;
    av_y+=y[i]/Np;
    av_z+=z[i]/Np;
  }  

  sprintf(filename,"coord_average_ex_num_%d.dat",ex_num);
  fp=fopen(filename,"a+");
  fprintf(fp,"%f\t%f\t%f\t%f\n",t,av_x,av_y,av_z);
  fclose(fp);

  return 0;

}
int all_particle_coord(double*x,double*y,double*z,double t,int Np,int ex_num ){
  int i;
  char filename[64];
  FILE *fp;
sprintf(filename,"coord_all_partricle_ex_num_%d.dat",ex_num);
  fp=fopen(filename,"a+");
  for(i=0;i<Np;i++){
    fprintf(fp,"%f\t%f\t%f\t%f\n",t,x[i],y[i],z[i]);
}
  fclose(fp);
  return 0;
}
int main()
{
  double t, avU=0.0,avK=0.0,Lx=0.0,Ly=0.0,Lz=0.0,Th;
  int i,Np,count=0;
  Np = 300;
  double r1=1.0,r2=1.4;
  double L;
  double* x, * y, * z,* vx, * vy,* vz,*a, * kx, * ky,* kz;
  x = new double[Np];
  y = new double[Np];
  z = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  vz = new double[Np];
  a = new double[Np];
  kx = new double[Np];
  ky = new double[Np];
  kz = new double[Np];
  double time_stable=1.e+2,time_eq=1.e+2;
  double sampling_time,time_stamp=0.;
  double sampling_time_max=1.e+3;
  double dt = 1.e-5, time_max = 1.e+3;
  char filename[128];
  int ex_num=10;
  L = pow (Np / 3.e-1 ,1.0/3.0);
  seed_change(ex_num);

  ini_box_rand(x, y, z,a, Np, L, r1, r2);//位置のboxsize依存を少しでも減らすために座標は-0.5*L~0.5*L
  ini(vx, vy,vz, Np);
  //HP
  Th = 0.0;
  //粒子をばらつかせる
  for (t = 0.; t < time_stable; t += dt) {
    
    calc_force_hs(x, y,z, L, Np, a, kx, ky,kz, &avU);
    eq_motion(x, y, z,vx, vy,vz, kx, ky,kz,dt, Np, &avK,&Lx,&Ly,&Lz, Th);
    p_bound(x, y,z, Np, L);
  }

  Th=1;//HPによる平衡化
 for (t = 0.; t < time_eq; t += dt) {
    calc_force_hs(x, y,z, L, Np, a, kx, ky,kz, &avU);
    eq_motion(x, y, z,vx, vy,vz, kx, ky,kz,dt, Np, &avK,&Lx,&Ly,&Lz, Th);
    p_bound(x, y,z, Np, L);
  }

//これ以降は相互作用なしなのでforceなし
ini_force(kx,ky,kz,Np);
//粒子の位置確認
//ini all_particle_coord_see(double *x,double *y,double *z,int Np,int ex_num){
 all_particle_coord_see(x,y,z,Np,ex_num,&avK,L);//位置のばらつきを見る。
//デバッグ用
printf("kx[1]は%f、ky[1]は%f、kz[1]は%f、0になるはず。",kx[1],ky[1],kz[1]);
sprintf(filename,"energy_time_B_20_ex_num_%d.txt",ex_num);
ofstream file;
file.open(filename);
double ave_x=0,ave_y=0,ave_z=0;
  sampling_time=5.*dt;
  int count2=0;
  //EoM
  for (t = 0; t < time_max; t += dt) {
    count++;
    count2++;
    non_equili_eq_motion(x, y, z,vx, vy, vz,kx,ky,kz,dt, Np, &avK,&Lx,&Ly,&Lz,L);
    //com_correction(x,y,z,&x_corr,&y_corr,&z_corr,Np,L);//独立なので重心とかない。
    p_bound(x, y, z,Np, L);
    if(count==1.e+6){
      output_M(&Lx,&Ly,&Lz,Np,t,ex_num,L);
      file <<t<<" "<< avK <<endl;
      output(x,y,z,vx,vy,vz,t,ex_num,Np);
      coord_average(x,y,z,Np,t,ex_num);
      count=0;
    }
    if(count2==5.e+7){
      all_particle_coord(x,y,z,t,Np,ex_num);
      count2=0;
}
    if(int(t/dt) == int((sampling_time + time_stamp)/dt)){
      output_t(x,y,z,Np,t,time_stamp,ex_num);
      sampling_time*=pow(10.,0.1);
      sampling_time=int(sampling_time/dt)*dt;
      if(sampling_time > sampling_time_max/pow(10.,0.1)){
	time_stamp=t;
	sampling_time=5.*dt;
      }
    }
  }
  file.close();
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] a;
  delete[] kx;
  delete[] ky;
  delete[] kz;
  return 0;
}
