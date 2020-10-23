# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
using namespace std;
# define Npm 300

double unif_rand(double left, double right)
{
  return left + (right - left) * rand() / RAND_MAX;
}

int ini_coord_rand(double* x, double* y,double* a, int Np, double L, double r1, double r2){
  int k, p;
  p = 0;
  for (k = 0; k < Np; k++) {
    x[k] = unif_rand(0, 1) * L;
    y[k] = unif_rand(0, 1) * L;
    //箱の中にランダムに配置
    if (p == 0) {
      a[k] = r1;
      p = p + 1;
      //異なる大きさの粒子を1対1でおいたプログラム、座標とセットでおいた
    }
    else {
      a[k] = r2;
      p = p - 1;
    }
  }
  return 0;
}


int ini(double* vx, double* vy,int Np) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
  }
  return 0;
}

//周期境界条件、tが更新されるたびに箱から出てないか確認
int p_bound(double* x, double* y,int Np, double L) {
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

  }
  return 0;
}

int calc_force_hs(double* x, double* y,double L, int Np, double* a, double* kx, double* ky,double* avU) {
  int i, j, k;
  *avU = 0.0;
  double r;
  double t, drU;
  double dx, dy;
  double aij;
  double cut;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
  }
  for (i = 0; i < Np; i++){
     for (j = i + 1; j < Np; j++){
       dx = x[i] - x[j];
       dy = y[i] - y[j];
      if (dx > (0.5 * L))
         dx -= L;
      if (dx < -(0.5 * L))
        dx += L;
      if (dy > (0.5 * L))
        dy -= L;
      if (dy < -(0.5 * L))
        dy += L ;
        aij = (a[i] + a[j]) / 2.0;
        r = sqrt(dx * dx + dy * dy);
        t = r / aij;
        cut = aij;
        if (r < cut) {
	  drU = -(1 - t) / aij; //analytical calculation of the 1'st derivative
        }
        else {
	  drU = 0.0;
	  continue;
        }
        kx[i] -= drU * dx / r; //reaction//
        kx[j] += drU * dx / r; //action//
        ky[i] -= drU * dy / r;
        ky[j] += drU * dy / r;
        *avU += (1 - t) * (1 - t);
      }
  }
  *avU /= double(Np);
  return 0;
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

int eq_motion(double* x, double* y,double* vx, double* vy, double dt, double* kx, double* ky,int Np, double* avK, double Th) {
  double zeta;
  zeta = 1.;
  int k;
  for (k = 0; k < Np; k++) {
    vx[k] += -vx[k] * zeta * dt + kx[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    vy[k] += -vy[k] * zeta * dt + ky[k] * dt + sqrt(2. * zeta * Th * dt) * gaussian_rand();
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k] ;
  }
  *avK = *avK / Np / 2.0;
  return 0;
}

int output(double *x,double *y,double *r1,int Np,double t,double dt){
  int i;
  char filename[128];//これはファイル名のための文字列//
  sprintf(filename,"bina2d_%.0f.dat",(t/dt));//ファイル名の決定//
  //%.0fは小数点0桁まで、%.1fは小数点一桁まで表示//
  ofstream file;
  file.open(filename);//ファイルを書き込み用にオープン//

   for(i=0;i<Np;i++)
    file << x[i] << " " << y[i]<< " "<< r1[i] << endl;
  file.close();
  return 0;
}


int main(){
  double t, avU=0.0, avK=0.0,Th;
  int i,count=0;
  double time_stable = 100.;//緩和時間
  double dt = 0.01;//刻み幅
  int Np = 300;
  double L;//箱
  double* x, * y, * vx, * vy,* a, * kx, * ky;//
  double r1=1.0, r2=1.4;//粒子の半径

  x = new double[Np];
  y = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  a = new double[Np];
  kx = new double[Np];
  ky = new double[Np];

  char filename[128];
  L = sqrt(double(Np) / 0.8);
  ini_coord_rand(x, y, a, Np, L, r1, r2);
  ini(vx, vy, Np);

  ofstream file;



  //HP
  Th = 0.0;

  for (t = 0.; t < time_stable; t += dt) {

  //緩和時間までdt刻みでUやKの時間平均をとったりしてる
    calc_force_hs(x, y, L, Np, a, kx, ky, &avU);
    eq_motion(x, y, vx, vy, dt, kx, ky, Np, &avK, Th);
    p_bound(x, y, Np, L);
  }

  for(t=dt;dt<time_stable;t +=dt) {
    count++;
    calc_force_hs(x, y, L, Np, a, kx, ky, &avU);
    eq_motion(x, y, vx, vy, dt, kx, ky, Np, &avK, Th);
    //avUとavKをadressにする意味について後で考える
    p_bound(x, y, Np, L);
    if(count==100){
      output(x,y,a,Np,t,dt);
      //これはいらんのでは、ファイル出力してないし、file <<t<<" "<< avU<<" "<< avK <<endl;
      count=0;
    }
  }

delete[] x;
delete[] y;
delete[] vx;
delete[] vy;
delete[] a;
delete[] kx;
delete[] ky;
return 0;
}
