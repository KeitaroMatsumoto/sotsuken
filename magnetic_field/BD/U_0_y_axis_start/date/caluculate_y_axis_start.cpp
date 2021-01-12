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
}

//初期条件を粒子が等間隔にy軸上に乗るようにした
int ini_box_rand(double* x, double* y, double* z,double* a, int Np, double L, double r1, double r2){
  int k, p;
  p = 0;
  double d = L/(Np+1);//粒子の配置する間隔
  for (k = 0; k < Np; k++) {
    x[k] = 0;
    y[k] = 0.5*L - d*(k+1);//上から順に配置
    z[k] = 0;//初期条件変更
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
    vx[j] = 600.0;//y軸にsetしてx軸方向に初速与える
    vy[j] = 0.0;
    vz[j] = 0.0;
  }
  return 0;
}




int eq_motion(double* x, double* y,double* z, double* vx, double* vy,double* vz,double dt,int Np, double* avK, double* px,double* py,double* pz,double* p_all,double Th) {
  double Z;
  double B;
  Z = 1.;
  B = 0;
  double Vx=0;
  double p=1;//特徴的な運動量スケール
  int k;
  *px = 0;
  *py = 0;
  *pz = 0;
  *p_all = 0;
  for (k = 0; k < Np; k++) {
    Vx = vx[k];//vyに代入するvxだけ新しいものになるのを防ぐため。
    vx[k] += -vx[k] * Z * dt  + vy[k] * B * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();//磁場による項を追加
    vy[k] += -vy[k] * Z * dt  -Vx * B * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();//vxだけ直後のt+dt秒後のやつになってるのを修正
    vz[k] += -vz[k] * Z * dt + sqrt(2. * Z * Th * dt) * gaussian_rand();
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
    if (x[k] < -0.5*L) {
      x[k] = x[k] + L;
    }
    if (x[k] > 0.5* L) {
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
      if (z[k] >  0.5*L) {
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
  double t,avK=0.0,Lx=0.0,Ly=0.0,Lz=0.0,Th;//軌道角運動量も追加
  double px=0.0,py=0.0,pz=0.0,p_all=0.0;//運動量も追加
  double x_corr=0.0,y_corr=0.0,z_corr=0.0;//重心のやつも追加
  double r1, r2;
  int i,Np,count=0;
  Np = 300;//今回は粒子の数が少ないからcell listなくても早く終えれた。
  r1 = 1.0; r2 = 1.4;
  double L;
  double* x, * y, * z,* vx, * vy,* vz, * a;
  x = new double[Np];
  y = new double[Np];
  z = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  vz = new double[Np];
  a = new double[Np];
  double dt = 0.0001, time_max = 50.; //parameters,dtは0.01
  double time_stable = 100.;//time_stableを30秒から100秒にしてみた
  char filename[128];
  int ex_num=16;
  seed_change(ex_num);
  L = pow (double(Np) / 0.8 ,1.0/3.0);
  ini_box_rand(x, y, z,a, Np, L, r1, r2);//ここからもう箱の座標中心を0にした
  ini(vx, vy,vz, Np);

  //HPによる粒子を散らばらせるプログラムもいらん。粒子の重なりはどうするか問題はある
  //BHHP
  Th = 1.0;
  int b=0;//ファイルに書く用、関数の引数にいれるの後で。
  double Vx=600;//ファイルに書き出す用のVxの初速

  sprintf(filename,"energy_time_B_%d_Np_%.0d_Th_%.0f_ini_v_%.0f_ex_num_%.0d_time_%.0f.txt",b,Np,Th,Vx,ex_num,time_max);
  ofstream file;
  file.open(filename);

  sprintf(filename,"orbital_magnetic_moment_B_%d_Np_%.0d_Th_%.0f_ini_v_%.0f_ex_num_%.0d_time_%.0f.txt",b,Np,Th,Vx,ex_num,time_max);//二つ目のファイル、軌道角運動量用に
  ofstream file2;
  file2.open(filename);

  sprintf(filename,"BD_center_of_gravity_B_%d_Np_%.0d_Th_%.0f_ini_v_%.0f_ex_nun_%.0d_time_%.0f.txt",b,Np,Th,Vx,ex_num,time_max);//３つ目のファイル、重心の位置を見る用に
  ofstream file3;
  file3.open(filename);

  sprintf(filename,"BD_momentum_B_%d_Np_%.0d_Th_%.0f_ini_v_%.0f_ex_num_%.0d_time_%.0f.txt",b,Np,Th,Vx,ex_num,time_max);//運動量の吐き出し
  ofstream file4;
  file4.open(filename);


  //EoM
  for (t = dt; t < time_max; t += dt) {
    count++;
    ///calc_force(x, y, z,L, Np, a, kx, ky, kz,&avU);
    eq_motion(x, y, z,vx, vy, vz,dt, Np, &avK,&px,&py,&pz,&p_all,Th);
    com_correction(x,y,z,&x_corr,&y_corr,&z_corr,Np,L);
    p_bound(x, y, z,Np, L);
    calc_M(x,y,z,vx,vy,vz,&Lx,&Ly,&Lz,L,Np);


    //(double* x, double* y,double* z,double* vx, double* vy, double* vz,double* Lx,double* Ly,double* Lz,double* del_Mx,double* del_My,double* del_Mz,double L)
    //int com_correction(double *x,double *y,double *z,double *x_corr,double *y_corr,double *z_corr,int Np,double L)
    if(count==1000){
      /// output(x,y,z,a,Np,t,dt);
      // output_v(vx,vy,vz,Np,t,dt);
      /// output_p(&px,&py,&pz,Np,t,dt);これももう吐き出さなくて良い、確認済み
      ///int output_p(double* px,double* py,double* pz,int Np,double t,double dt)
      /// output_f(kx,ky,kz,Np,t,dt);

      //int output_v(double *vx,double *vy,double *vz)



      file <<t<<" "<< avK <<endl;
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
  return 0;
}
