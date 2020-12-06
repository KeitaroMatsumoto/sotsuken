# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
# define Npm 300
using namespace std;
int output(double *x,double *y,double *z,double *r1,int Np,double t,double dt){
  int i;
  char filename[128];
  sprintf(filename,"bina2d_%.0f.dat",(t/dt));
  ofstream file;
  file.open(filename);
  for(i=0;i<Np;i++)
    file << x[i] << " " << y[i]<< " " <<z[i]<<" "<< r1[i] << endl;
  file.close();
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

int ini_box_rand(double* x, double* y, double* z,double* a, int Np, double L, double r1, double r2){
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

int ini(double* vx, double* vy, double* vz,int Np) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
    vz[j] = 0.0;
  }
  return 0;
}

int calc_force_hs(double* x, double* y,double* z, double L, int Np, double* a, double* kx, double* ky, double* kz,double* avU) {
  int i, j, k;
  *avU = 0.0;
  double r;
  double t, f;
  double dx, dy,dz;
  double aij;
  double cut;
  for (k = 0; k < Np; k++){
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


int calc_force(double* x, double* y, double* z,double L, int Np, double* a, double* kx, double* ky,double* kz, double* avU){ //力の計算
  int i, j, k;
  *avU = 0.0;
  double r;
  double w,w2,w4,w12,f;
  double dx, dy,dz;
  double aij;
  double cut;
  cut = 3.0;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
    kz[k] = 0.0;
  }
  for (j = 0; j < Np; j++)
    {
      for (i = j + 1; i < Np; i++)
	{
          dx = x[j] - x[i];
          dy = y[j] - y[i];
          dz = z[j] - z[i];
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
          aij = (a[j] + a[i]) / 2.0;
          r = sqrt(dx * dx + dy * dy + dz * dz);
          w = aij / r;
	  w2=w*w;
	  w4=w2*w2;
	  w12=w4*w4*w4;
          if (r < cut) {
            f = (-12.0) * w12 / r; //analytical calculation of the 1'st derivative
          }
          else {
            f = 0.0;
            continue;
          }
          kx[j] -= f * dx / r;
          kx[i] += f * dx / r;
          ky[j] -= f * dy / r;
          ky[i] += f * dy / r;
          kz[j] -= f * dz / r;
          kz[i] += f * dz / r;
          *avU += w12;
        }
    }
    *avU /= double(Np);
    return 0;
  }


int eq_motion(double* x, double* y,double* z, double* vx, double* vy,double* vz,double dt, double* kx, double* ky, double* kz,int Np, double* avK, double* px,double* py,double* pz,double* p_all,double Th) {
  double Z;
  double B;
  Z = 1.;
  B = 2;
  double Vx=0;
  double p=1;//特徴的な運動量スケール
  int k;
  *px = 0;
  *py = 0;
  *pz = 0;
  *p_all = 0;
  for (k = 0; k < Np; k++) {
    Vx = vx[k];//vyに代入するvxだけ新しいものになるのを防ぐため。
    vx[k] += -vx[k] * Z * dt + kx[k] * dt + vy[k] * B * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();//磁場による項を追加
    vy[k] += -vy[k] * Z * dt + ky[k] * dt -Vx * B * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();//vxだけ直後のt+dt秒後のやつになってるのを修正
    vz[k] += -vz[k] * Z * dt + kz[k] * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();
    ///無次元化された質量は1である
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    z[k] += vz[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k];
    *px += p * vx[k];
    *py += p * vy[k];
    *pz += p * vz[k];//個別の運動量p*=v*
    *p_all += vx[k] + vy[k] + vz[k];
  }
  *avK = *avK / Np / 2.0;
  return 0;
}

int calc_M(double* x, double* y,double* z,double* vx, double* vy, double* vz,double* Lx,double* Ly,double* Lz,double L,int Np){
  double μ;
  μ = 1.;//無次元化した際に出てくる係数μ*

  int k;
  * Lx=0;
  *Ly=0;
  *Lz=0;
  for(k=0;k<Np;k++){
    *Lx +=  μ * (y[k] * vz[k] - z[k] * vy[k]);
    *Ly +=  μ * (z[k] * vx[k] - x[k] * vz[k]);
    *Lz +=  μ * (x[k] * vy[k] - y[k] * vx[k]);
  }
  return 0;
}

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
      if (z[k] >  L) {
        z[k] = z[k] - L;
    }
  }
  return 0;
}

//重心を記録するプログラムを追加
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

    *x_corr+=dx/Np; //center of mass displacement.x,重心がどれほど動いたかを記録
    *y_corr+=dy/Np;
    *z_corr+=dz/Np;


    x0[i]=x[i];
    y0[i]=y[i];
    z0[i]=z[i];
  }
  return 0;
}

int main(){
  double t, avU=0.0, avK=0.0,Lx=0.0,Ly=0.0,Lz=0.0,Th;//軌道角運動量も追加
  double del_Mx=0.0,del_My=0.0,del_Mz = 0.0;//Mゆらぎ
  double px=0.0,py=0.0,pz=0.0,p_all=0.0;//運動量も追加
  double x_corr=0.0,y_corr=0.0,z_corr=0.0;//重心のやつも追加
  double r1, r2;
  int i,Np,count=0;
  Np = 300;//今回は粒子の数が少ないからcell listなくても早く終えれた。
  r1 = 1.0; r2 = 1.4;
  double L;
  double* x, * y, * z,* vx, * vy,* vz, * a, * kx, * ky,* kz;
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
  double dt = 0.01, time_max = 100.; //parameters,dtは0.01
  double time_stable = 100.;//time_stableを30秒から100秒にしてみた
  char filename[128];
  L = pow (double(Np) / 0.8 ,1.0/3.0);
  ini_box_rand(x, y, z,a, Np, L, r1, r2);
  ini(vx, vy,vz, Np);

  //HP
  Th = 0.0;

  for (t = 0.; t < time_stable; t += dt) {
    calc_force_hs(x,y,z,L,Np,a, kx, ky,kz, &avU);
    eq_motion(x, y, z,vx, vy,vz, dt, kx, ky,kz, Np, &avK,&px,&py,&pz,&p_all, Th);
    p_bound(x, y,z, Np, L);
  }

  //BHHP
  Th = 1.0;

  sprintf(filename,"energy_time2_stabletime_B_2_100.txt");
  ofstream file;
  file.open(filename);

  sprintf(filename,"orbital_magnetic_moment_B_0_stabletime_100.txt");//二つ目のファイル、軌道角運動量用に
  ofstream file2;
  file2.open(filename);

  sprintf(filename,"BD_center_of_gravity_B_0.txt");//３つ目のファイル、重心の位置を見る用に
  ofstream file3;
  file3.open(filename);

  sprintf(filename,"BD_momentum_B_0.txt");//運動量の吐き出し
  ofstream file4;
  file4.open(filename);

  //EoM
  for (t = dt; t < time_max; t += dt) {
    count++;
    calc_force(x,y,z,L,Np,a,kx, ky, kz,&avU);
    eq_motion(x, y, z,vx, vy, vz,dt, kx, ky,kz, Np, &avK,&px,&py,&pz,&p_all,Th);
    com_correction(x,y,z,&x_corr,&y_corr,&z_corr,Np,L);
    p_bound(x, y, z,Np, L);
    calc_M(x,y,z,vx,vy,vz,&Lx,&Ly,&Lz,L,Np);

    //remove_boxsize_bound(double* x,double* y,double* z,double* x_l,double* y_l,double* z_l,int Np,double L){
    ///int calc_del_M(double* x_l,double* y_l,double* z_l,double* vx,double* vy,double* vz,double* del_Mx,double* del_My,double* del_Mz,int Np
    ///calc_Mではdel_Mxなどはやらない、除く
    ///remove_bound___をしてから、x_lなどを用いてdel_Mなどをけいさんする関数をいれるcalc_del_M

    //(double* x, double* y,double* z,double* vx, double* vy, double* vz,double* Lx,double* Ly,double* Lz,double* del_Mx,double* del_My,double* del_Mz,double L)
    //int com_correction(double *x,double *y,double *z,double *x_corr,double *y_corr,double *z_corr,int Np,double L)
    if(count==100){
      /// output(x,y,z,a,Np,t,dt);
      file <<t<<" "<< avU<<" "<< avK <<" "<<avU+avK<<endl;
      file2<<t<<" "<<Lx<<" "<<Ly<<" "<<Lz<<" "<<Lx+Ly+Lz<<endl;
      file3<<t<<" "<<x_corr<<" "<<y_corr<<" "<<z_corr<<" "<<endl;
      file4<<t<<" "<<px<<" "<<py<<" "<<pz<<" "<<p_all<<" "<<endl;
      count=0;
    }

  }
  file.close();
  file2.close();
  file3.close();
  file4.close();

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
