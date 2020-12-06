# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <fstream>
# include <time.h>
using namespace std;
# define Npm 300
int main(){
  int i,j,k;
  int Np = 300;
  double* x, * y, * z, * r1;
  x = new double[Np];
  y = new double[Np];
  z = new double[Np];
  r1 = new double[Np];

  char filename[64];
  ifstream file;//ファイルの読み取り

  sprintf(filename,"bina2d_9900.dat");
  file.open(filename);

  for(i=0;i<Np;i++){
    file >> x[i] >> y[i] >> z[i] >> r1[i];
  }
  file.close();

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] r1;
  return 0;
}
