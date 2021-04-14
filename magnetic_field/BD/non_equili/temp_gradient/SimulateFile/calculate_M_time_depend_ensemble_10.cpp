#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
int main(){
char filename1[200];
ifstream file1;//読み取り
char filename2[200];
ifstream file2;//読み取り
char filename3[200];
ifstream file3;//読み取り
char filename4[200];
ifstream file4;//読み取り
char filename5[200];
ifstream file5;//読み取り
char filename6[200];
ifstream file6;//読み取り
char filename7[200];
ifstream file7;//読み取り
char filename8[200];
ifstream file8;//読み取り
char filename9[200];
ifstream file9;//読み取り
char filename10[200];
ifstream file10;//読み取り


int ex_num=10;//読み取るアンサンブル数
double max=1.e+2;//データの要素数


char filename11[200];
ofstream file11;//書き出し
char filename12[200];
ofstream file12;//書き出し

double Vx_1[10000],Vy_1[10000],Vz_1[10000],t_1[10000];
double Vx_2[10000],Vy_2[10000],Vz_2[10000],t_2[10000];
double Vx_3[10000],Vy_3[10000],Vz_3[10000],t_3[10000];
double Vx_4[10000],Vy_4[10000],Vz_4[10000],t_4[10000];
double Vx_5[10000],Vy_5[10000],Vz_5[10000],t_5[10000];//読み取ったデータ用,データの要素数に合わせたほうがいいかも
double Vx_6[10000],Vy_6[10000],Vz_6[10000],t_6[10000];
double Vx_7[10000],Vy_7[10000],Vz_7[10000],t_7[10000];
double Vx_8[10000],Vy_8[10000],Vz_8[10000],t_8[10000];
double Vx_9[10000],Vy_9[10000],Vz_9[10000],t_9[10000];
double Vx_10[10000],Vy_10[10000],Vz_10[10000],t_10[10000];

double ave_Vx[10000],ave_Vy[10000],ave_Vz[10000],t[10000];//吐き出し用
double sigma_Vx[10000],sigma_Vy[10000],sigma_Vz[10000];

int i;

double time_eq=0;

for(i=0;i<max;i++){
  ave_Vx[i]=0.;
  ave_Vy[i]=0.;
  ave_Vz[i]=0.;
  sigma_Vx[i]=0.;
  sigma_Vy[i]=0.;
  sigma_Vz[i]=0.;//初期化
}

sprintf(filename1,"orbital_magnetic_moment_Np_300_ex_num_1.dat");
file1.open(filename1);
for(i=0;i<max;i++){
  file1 >>t_1[i] >>Vx_1[i]>>Vy_1[i]>>Vz_1[i];
}
file1.close();

sprintf(filename2,"orbital_magnetic_moment_Np_300_ex_num_2.dat");
file2.open(filename2);
for(i=0;i<max;i++){
  file2 >>t_2[i] >>Vx_2[i]>>Vy_2[i]>>Vz_2[i];
}
file2.close();

sprintf(filename3,"orbital_magnetic_moment_Np_300_ex_num_3.dat");
file3.open(filename3);
for(i=0;i<max;i++){
  file3 >>t_3[i] >>Vx_3[i]>>Vy_3[i]>>Vz_3[i];
}
file3.close();

sprintf(filename4,"orbital_magnetic_moment_Np_300_ex_num_4.dat");
file4.open(filename4);
for(i=0;i<max;i++){
  file4 >>t_4[i] >>Vx_4[i]>>Vy_4[i]>>Vz_4[i];
}
file4.close();

sprintf(filename5,"orbital_magnetic_moment_Np_300_ex_num_5.dat");
file5.open(filename5);
for(i=0;i<max;i++){
  file5 >>t_5[i] >>Vx_5[i]>>Vy_5[i]>>Vz_5[i];
}
file5.close();
sprintf(filename6,"orbital_magnetic_moment_Np_300_ex_num_6.dat");
file6.open(filename6);
for(i=0;i<max;i++){
  file6 >>t_6[i] >>Vx_6[i]>>Vy_6[i]>>Vz_6[i];
}
file6.close();

sprintf(filename7,"orbital_magnetic_moment_Np_300_ex_num_7.dat");
file7.open(filename7);
for(i=0;i<max;i++){
  file7 >>t_7[i] >>Vx_7[i]>>Vy_7[i]>>Vz_7[i];
}
file7.close();

sprintf(filename8,"orbital_magnetic_moment_Np_300_ex_num_8.dat");
file8.open(filename8);
for(i=0;i<max;i++){
  file8 >>t_8[i] >>Vx_8[i]>>Vy_8[i]>>Vz_8[i];
}
file8.close();

sprintf(filename9,"orbital_magnetic_moment_Np_300_ex_num_9.dat");
file9.open(filename9);
for(i=0;i<max;i++){
  file9 >>t_9[i] >>Vx_9[i]>>Vy_9[i]>>Vz_9[i];
}
file9.close();

sprintf(filename10,"orbital_magnetic_moment_Np_300_ex_num_10.dat");
file10.open(filename10);
for(i=0;i<max;i++){
  file10 >>t_10[i] >>Vx_10[i]>>Vy_10[i]>>Vz_10[i];
}
file10.close();
//デバグ用
 for(i=0;i<max;i++){
   printf("%f,%f, %f %f,%f,%f,%f,%f,%f,%f\n",Vz_1[i],Vz_2[i],Vz_3[i],Vz_4[i],Vz_5[i],Vz_6[i],Vz_7[i],Vz_8[i],Vz_9[i],Vz_10[i]);
 }

//各時間の物理量の平均
for(i=0;i<max;i++){
  ave_Vx[i]=(Vx_1[i]+Vx_2[i]+Vx_3[i]+Vx_4[i]+Vx_5[i]+Vx_6[i]+Vx_7[i]+Vx_8[i]+Vx_9[i]+Vx_10[i])/(ex_num);
  ave_Vy[i]=(Vy_1[i]+Vy_2[i]+Vy_3[i]+Vy_4[i]+Vy_5[i]+Vy_6[i]+Vy_7[i]+Vy_8[i]+Vy_9[i]+Vy_10[i])/(ex_num);
  ave_Vz[i]=(Vz_1[i]+Vz_2[i]+Vz_3[i]+Vz_4[i]+Vz_5[i]+Vz_6[i]+Vz_7[i]+Vz_8[i]+Vz_9[i]+Vz_10[i])/(ex_num);

}

for(i=time_eq;i<max;i++){
  sigma_Vx[i]+= (Vx_1[i]-ave_Vx[i])*(Vx_1[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_2[i]-ave_Vx[i])*(Vx_2[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_3[i]-ave_Vx[i])*(Vx_3[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_4[i]-ave_Vx[i])*(Vx_4[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_5[i]-ave_Vx[i])*(Vx_5[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_6[i]-ave_Vx[i])*(Vx_6[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_7[i]-ave_Vx[i])*(Vx_7[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_8[i]-ave_Vx[i])*(Vx_8[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_9[i]-ave_Vx[i])*(Vx_9[i]-ave_Vx[i])/(ex_num);
  sigma_Vx[i]+= (Vx_10[i]-ave_Vx[i])*(Vx_10[i]-ave_Vx[i])/(ex_num);

  sigma_Vy[i]+= (Vy_1[i]-ave_Vy[i])*(Vy_1[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_2[i]-ave_Vy[i])*(Vy_2[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_3[i]-ave_Vy[i])*(Vy_3[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_4[i]-ave_Vy[i])*(Vy_4[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_5[i]-ave_Vy[i])*(Vy_5[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_6[i]-ave_Vy[i])*(Vy_6[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_7[i]-ave_Vy[i])*(Vy_7[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_8[i]-ave_Vy[i])*(Vy_8[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_9[i]-ave_Vy[i])*(Vy_9[i]-ave_Vy[i])/(ex_num);
  sigma_Vy[i]+= (Vy_10[i]-ave_Vy[i])*(Vy_10[i]-ave_Vy[i])/(ex_num);

  sigma_Vz[i]+= (Vz_1[i]-ave_Vz[i])*(Vz_1[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_2[i]-ave_Vz[i])*(Vz_2[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_3[i]-ave_Vz[i])*(Vz_3[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_4[i]-ave_Vz[i])*(Vz_4[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_5[i]-ave_Vz[i])*(Vz_5[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_6[i]-ave_Vz[i])*(Vz_6[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_7[i]-ave_Vz[i])*(Vz_7[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_8[i]-ave_Vz[i])*(Vz_8[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_9[i]-ave_Vz[i])*(Vz_9[i]-ave_Vz[i])/(ex_num);
  sigma_Vz[i]+= (Vz_10[i]-ave_Vz[i])*(Vz_10[i]-ave_Vz[i])/(ex_num);

  sigma_Vx[i] = pow(sigma_Vx[i],0.5);
  sigma_Vy[i] = pow(sigma_Vy[i],0.5);
  sigma_Vz[i] = pow(sigma_Vz[i],0.5);

  sigma_Vx[i]=sigma_Vx[i]/pow(10,0.5);
  sigma_Vy[i]=sigma_Vy[i]/pow(10,0.5);
  sigma_Vz[i]=sigma_Vz[i]/pow(10,0.5);//アンサンブル数によるゆらぎの減少
}
//デバグ用
 for(i=0;i<max;i++){
   printf("t_1[%d]=%f、sigma_Vx[%d]=%f、sigma_Vy[%d]=%f、sigma_Vz[%d]=%f \n",i,t_1[i],i,sigma_Vx[i],i,sigma_Vy[i],i,sigma_Vz[i]);
 }


sprintf(filename11,"M_ave_and_variance_B_3_ensemble_num_10.dat");
file11.open(filename11);
sprintf(filename12,"M_real_error_var_value_date_ensemble_num_10.dat");
file12.open(filename12);
for(i=0;i<max;i++){
  file11<<t_1[i]<<" "<<ave_Vx[i]<<" "<<ave_Vy[i]<<" "<<ave_Vz[i]<<" "<<sigma_Vx[i]<<" "<<sigma_Vy[i]<<" "<<sigma_Vz[i]<<endl;
  file12<<t_1[i]<<" "<<-sigma_Vx[i]<<" "<<sigma_Vx[i]<<" "<<-sigma_Vy[i]<<" "<<sigma_Vy[i]<<" "<<-sigma_Vz[i]<<" "<<sigma_Vz[i]<<endl;
}
file11.close();
file12.close();

  return 0;
}
