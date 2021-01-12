
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
int main(){
char filename[200];
ifstream file;//読み取り
char filename2[200];
ofstream file2;//書き出し
int ex_num=5;//読み取るアンサンブルのインデント
double Mx[10000],My[10000],Mz[10000],t[10000];
double ave_Mx=0,ave_My=0,ave_Mz=0;
double sigma_Mx=0,sigma_My=0,sigma_Mz=0;
int i;
double max=1.e+4;//今回dt=1.e-3でcount1000で吐き出し、timemaxが１.e+4なので要素数は10000でいい
double time_eq=10;//平衡化してから計算

sprintf(filename,"orbital_magnetic_moment_B_0_Np_300_Th_1_ini_v_600_ex_num_%d_time_10000.txt",ex_num);
file.open(filename);
for(i=0;i<max;i++){
  file >>t[i] >>Mx[i]>>My[i]>>Mz[i];
}
file.close();


for(i=time_eq;i<max;i++){
  ave_Mx+=Mx[i]/(max-time_eq);
  ave_My+=My[i]/(max-time_eq);
  ave_Mz+=Mz[i]/(max-time_eq);

}

for(i=time_eq;i<max;i++){
  sigma_Mx+= (Mx[i]-ave_Mx)*(Mx[i]-ave_Mx)/(max-time_eq);
  sigma_My+= (My[i]-ave_My)*(My[i]-ave_My)/(max-time_eq);
  sigma_Mz+= (Mz[i]-ave_Mz)*(Mz[i]-ave_Mz)/(max-time_eq);

}
sigma_Mx = pow(sigma_Mx,0.5);
sigma_My = pow(sigma_My,0.5);
sigma_Mz = pow(sigma_Mz,0.5);

sprintf(filename2,"M_ave_and_variance_ex_num_%d.dat",ex_num);
file2.open(filename2);
file2<<ave_Mx<<" "<<ave_My<<" "<<ave_Mz<<" "<<sigma_Mx<<" "<<sigma_My<<" "<<sigma_Mz<<endl;
file2.close();

  return 0;
}
