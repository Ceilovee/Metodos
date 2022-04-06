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
// m+1 radios y n angulos
// espero que el vector te tenga n valores
void altoHorno(int ri, int re, int m, int n, double isoterma, double ti, vector<double> te){
    // tengo que:
    // A= 1 0 0 0  segun se mjultiplique x te
    // b= (ti(x n), 0 ... 0, te)
    // x= (tm=1n=1, ....... tm=m+1n=n)
    vector<double> b(m*n);
    for(int i=0; i<b.size(); i++){
        if(i<n){
            b[i]=ti;
        }else if(i>=(b.size()-n)){
            b[i]=te[n-(b.size()-i)];
        }else{
            b[i]=0;
        }
    }
    // con esto se crea b de la forma anterior

    // creo la matriz con los coeficientes
    vector<vector<double>> A;
    int i=0;
    vector<double>aux(n*m,0);
    double c1 = (1/(((re-ri)/m)^2) - (1/((re-ri)*((re-ri)/m)))); // ????
    double c2 = (-2/(((re-ri)/m)^2))+(1/((re-ri)*((re-ri)/m)))-(2/((1))); // esto esta mal
    double c3 = (1/(((re-ri)/m)^2));
    double c4 = 0;
    double c5 = 0;
    for(int radio=0; radio<m; radio++){
        for(int angulo=0; angulo<n; angulo++){
            // esto se hace por cada fila de la matriz
            aux[i]=c2;
            // por como esta ordenado el vector b y x,
            // la matriz es cada fila una temperatura con angulo k y radio r
            // donde la fila anterior es la temp de angulo k-1 y radio r
            // si k<0 entonces es la temp con radio r-1
            // de la misma forma quedan organizadas las temp dentro de las filas
            // t(k,r-1) .... t(k-1,r) t(k,r) t(k+1,r) .... t(k,r+1) ....
            A[i]=aux;
            aux.clear();
            i++;
        }
    }

    // despues deberiamos calcular la x con las factorizaciones anteriores
    vector<double> x= resolverTriangularxLU(A,b);

}

int main() {
    std::cout << "Hello, World!" << std::endl;
    //vector<vector<double>> v = {{2,1,-1,3},{-2,0,0,0},{4,1,-2,4},{-6,-1,2,-3}};
   // mostrarMatriz(facLU(v));
   // vector<double> b={13,-2,24,-10};
   // vector<double> x= (resolverTriangularxLU(facLU(v),b));
    altoHorno(10,10,5,5,50,100,{20,15,18,16,25});
    return 0;
}

