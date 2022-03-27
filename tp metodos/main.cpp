#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <stdlib.h>

using namespace std;

void mostrarMatriz(vector<vector<double>> m){
    for(int i=0;i<m.size();i++){
        for(int j=0;j<m.size();j++) {
            cout<< m[i][j]<<" ";
        }
        cout<< endl;
    }
}

vector<vector<double>> crearRandom(int n){
    vector<vector<double>> res(n, vector<double>(n));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){ //  despues poner rango de 1499-0
            //double rand= rand(); //hacer numero random
            //res[i][j]= rand;
        }
    }
    return res;
}

vector<vector<double>> eliminacionGauss(vector<vector<double>> m){
    for(int i=0;i<m.size();i++){
        // if es 0 mover fila
        for(int j=i+1;j<m.size();j++){
            double num =m[j][i]/m[i][i];
            for(int k=i;k<m.size()+1;k++){
                m[j][k] = m[j][k] - (num* m[i][k]);
            }
        }
    }
    return m;
}

vector<vector<double>> facLU(vector<vector<double>> m){
    for(int i=0;i<m.size();i++){
        // if es 0 mover fila
        for(int j=i+1;j<m.size();j++){
            double num =m[j][i]/m[i][i];
            m[j][i]=num;
            for(int k=i+1;k<m.size()+1;k++){
                m[j][k] = m[j][k] - (num* m[i][k]);
            }
        }
    }
    return m;
}

vector<double> multiplicarL(vector<vector<double>> m, vector<double> b){
    // esto no funciona ARREGLAR
    vector<double> aux=b;
    for(int i=0;i<m.size();i++){
        for(int j=0;j<i;j++){
            b[j]+=m[i][j]*b[i];
        }
    }
    return b;
}

vector<double> resolverTriangularxLU(vector<vector<double>> m,vector<double> b){
    vector<double> x(b.size());
    b= multiplicarL(m,b);
    for(int i= m.size()-1; 0<=i; i--){
        for(int j=i+1; j<m.size();j++){
            b[i]-= m[i][j]*x[j];
        }
        x[i]=b[i]/m[i][i];
    }
    return x;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    vector<vector<double>> v = {{2,1,-1,3},{-2,0,0,0},{4,1,-2,4},{-6,-1,2,-3}};
    mostrarMatriz(facLU(v));
    vector<double> b={13,-2,24,-10};
    vector<double> x= (resolverTriangularxLU(facLU(v),b));
    return 0;
}

