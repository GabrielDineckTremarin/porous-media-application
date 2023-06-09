#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iomanip>
using namespace std;

//https://pt.stackoverflow.com/questions/108998/warning-em-compara%c3%a7%c3%a3o-entre-floats-como-proceder
bool isEqual(double a, double b) {
//cout << fixed << setprecision(8);
    return abs(a - b) <= 1e-5 * abs(a);
}

int main(int argc, char** argv){

  string arq1="saida200_OK.dat";
  //string arq1="teste200.dat";
  //string arq1="saida_200_2000.dat";
  //string arq1="saida51.dat";
  string arq2="data/output_variables.dat";
  ifstream arq(arq1);
  ifstream arqF(arq2);

  if((!arq) || (!arqF)){
    cout << "Não foi possivel abrir o arquivo" << endl;
    return 0;
  }

  int imax=51, erro=0;
  int jmax=63;
  double X[imax][jmax], Y[imax][jmax], u[imax][jmax], v[imax][jmax],T[imax][jmax];
  double Z[imax][jmax], H[imax][jmax], Yi[imax][jmax];
  double Xb[imax][jmax], Yb[imax][jmax], ub[imax][jmax], vb[imax][jmax],Tb[imax][jmax];
  double Zb[imax][jmax], Hb[imax][jmax], Yib[imax][jmax];
  cout << fixed << setprecision(15);

  for(int i=0;i<imax;i++){
	for(int j=0;j<jmax;j++){
          arq >> X[i][j]; //cout << X[i][j] << endl;
          arq >> Y[i][j];
          arq >> u[i][j];
          arq >> v[i][j];
          arq >> T[i][j];
          arq >> Z[i][j];
          arq >> H[i][j];
          arq >> Yi[i][j];
	}
  }

  arq.close();
  for(int i=0;i<imax;i++){
        for(int j=0;j<jmax;j++){
          arqF >> Xb[i][j];
          arqF >> Yb[i][j];
          arqF >> ub[i][j];
          arqF >> vb[i][j];
          arqF >> Tb[i][j];
          arqF >> Zb[i][j];
          arqF >> Hb[i][j];
          arqF >> Yib[i][j];
        }
  }  
  arqF.close();

  for(int i=0;i <imax;i++){
        for(int j=0; j<jmax;j++){
          if(!isEqual(Xb[i][j], X[i][j])){ cout << "Error X "<< Xb[i][j] << " != " << X[i][j] << endl; erro+=1;}
          if(!isEqual(Yb[i][j], Y[i][j])){ cout << "Error Y "<< Yb[i][j] << " != " << Y[i][j] << endl; erro+=1;}
          if(!isEqual(ub[i][j], u[i][j])){ cout << "Error U "<< ub[i][j] << " != " << u[i][j] << endl; erro+=1;}
          if(!isEqual(vb[i][j], v[i][j])){ cout << "Error V "<< vb[i][j] << " != " << v[i][j] << endl; erro+=1;}
          if(!isEqual(Tb[i][j], T[i][j])){ cout << "Error T "<< Tb[i][j] << " != " << T[i][j] << endl; erro+=1;}
          if(!isEqual(Zb[i][j], Z[i][j])){ cout << "Error Z "<< Zb[i][j] << " != " << Z[i][j] << endl; erro+=1;}
          if(!isEqual(Hb[i][j], H[i][j])){ cout << "Error H "<< Hb[i][j] << " != " << H[i][j] << endl; erro+=1;}
          if(!isEqual(Yib[i][j],Yi[i][j])){ cout << "Error Yi "<< Yib[i][j] << " != " <<Yi[i][j] << endl; erro+=1;}
        }
  }  

if(erro==0){cout << "Arquivos " << arq1 << " e " << arq2 << " são iguais" << endl;}

  return 0;
}
