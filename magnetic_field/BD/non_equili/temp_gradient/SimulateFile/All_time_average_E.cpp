#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
int main(){
char filename1[200];
ifstream file1;//書き出し


double max=1.e+2;//データの要素数


char filename6[200];
ofstream file6;//書き出し

double Vx_1[10000],Vy_1[10000],Vz_1[10000],t_1[10000];
double sigma_x_1[10000],sigma_y_1[10000],sigma_z_1[10000];
double ave_Vx=0,ave_Vy=0,ave_Vz=0;//吐き出し用

int i;


sprintf(filename1,"M_ave_and_variance_B_3_ensemble_num_10.dat");
file1.open(filename1);
for(i=0;i<max;i++){
  file1 >>t_1[i] >>Vx_1[i]>>Vy_1[i]>>Vz_1[i]>>sigma_x_1[i]>>sigma_y_1[i]>>sigma_z_1[i];
}
file1.close();


//デバグ用
 printf("t_1[%d]=%f、Vx_1[%d]=%f \n",0,t_1[0],0,Vx_1[0]);


//各時間の物理量の平均
for(i=0;i<max;i++){
  ave_Vx+=Vx_1[i]/(max);
  ave_Vy+=Vy_1[i]/(max);
  ave_Vz+=Vz_1[i]/(max);

}





sprintf(filename6,"All_time_orbital_magnetic_momentum_B_0_ensemble_10.dat");
file6.open(filename6);
file6<<ave_Vx<<" "<<ave_Vy<<" "<<ave_Vz<<endl;
file6.close();

  return 0;
}
