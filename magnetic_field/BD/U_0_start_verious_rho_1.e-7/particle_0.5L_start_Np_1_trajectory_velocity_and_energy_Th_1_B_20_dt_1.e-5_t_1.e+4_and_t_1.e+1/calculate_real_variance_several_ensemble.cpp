#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
int main(){
  char filename[200];
  ifstream file;
  char filename2[200];
  ifstream file2;
  char filename3[200];
  ifstream file3;
  char filename4[200];
  ifstream file4;
  char filename5[200];
  ifstream file5;
  char filename6[200];
  ofstream file6;
  int ex_num;
  double K1[10000],K2[10000],K3[10000],K4[10000],K5[10000],t1[10000],t2[10000],t3[10000],t4[10000],t5[10000];
  double av_K=0;
  double sigma_K=0;
  int i,j;
  double time_max=1.e+4;
  double time_eq=10;

  sprintf(filename,"energy_time_B_0_ex_num_1.txt");
  file.open(filename);
   for(i=0;i<time_max;i++){
      file>>t1[i]>>K1[i];
}
    file.close();

 sprintf(filename2,"energy_time_B_0_ex_num_2.txt");
 file2.open(filename2);
   for(i=0;i<time_max;i++){
      file2>>t2[i]>>K2[i];
}
    file2.close();

 sprintf(filename3,"energy_time_B_0_ex_num_3.txt");
 file3.open(filename3);
   for(i=0;i<time_max;i++){
      file3>>t3[i]>>K3[i];
}
    file3.close();

 sprintf(filename4,"energy_time_B_0_ex_num_4.txt");
 file4.open(filename4);
   for(i=0;i<time_max;i++){
      file4>>t4[i]>>K4[i];
}
    file4.close();

 sprintf(filename5,"energy_time_B_0_ex_num_5.txt");
 file5.open(filename5);
   for(i=0;i<time_max;i++){
      file5>>t5[i]>>K5[i];
}
    file5.close();


    for(i=time_eq;i<time_max;i++){
      av_K+=(K1[i]+K2[i]+K3[i]+K4[i]+K5[i])/(time_max - time_eq)/5;
}

    for(i=time_eq;i<time_max;i++){
      sigma_K+=(K1[i]-av_K)*(K1[i]-av_K)+(K2[i]-av_K)*(K2[i]-av_K)+(K3[i]-av_K)*(K3[i]-av_K)+(K4[i]-av_K)*(K4[i]-av_K)+(K5[i]-av_K)*(K5[i]-av_K);
}
    sigma_K=sigma_K/(time_max-time_eq)/5;//これは分散
    sigma_K=pow(sigma_K,0.5);//これはゆらぎ

    sprintf(filename6,"av_K_and_variance_ex_num_from_1_to_5.dat");
    file6.open(filename6);
    file6<<av_K<<" "<<sigma_K<<endl;
    file6.close();
    printf("K1は%f,K2は%f,K3は%f、K4は%f,K5は%f",K1[2],K2[2],K3[2],K4[2],K5[2]);

return 0;
}
