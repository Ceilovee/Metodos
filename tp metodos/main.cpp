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

/*
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
*/
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


vector<double> resolverTriangular_L(vector<vector<double>> m, vector<double> b){
    int n= m.size();
    vector<double> x(n);  
    
    for (int i = 0; i < n; i++){
        x[i] = b[i];
        
        for (int j=0; j<i; j++){
            x[i] -= m[i][j]*x[j];
        }
    }
 
    printf("\nSolution for the system:\n");
    for (int i=0; i<n; i++)
        printf("%lf\n", x[i]);
    
    return x;
}
 
vector<double> resolverTriangular(vector<vector<double>> m, vector<double> b){
    int n= m.size();
    vector<double> x(n);  
    
    for (int i = n-1; i >= 0; i--){
        x[i] = b[i];
 
        for (int j=i+1; j<n; j++){
            x[i] -= m[i][j]*x[j];
        }
        x[i] = x[i]/m[i][i];
    }
 
    printf("\nSolution for the system:\n");
    for (int i=0; i<n; i++)
        printf("%lf\n", x[i]);
    
    return x;
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
    b = resolverTriangular_L(m,b);
    x = resolverTriangular(m,b);
    return x;
}

vector<double> armarB(int n, int m, vector<double> te, double ti ){
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
    return b;
}

void encontrarIsoterma(vector<double> x, double isoterma){
    for(int i=0; i < x.size();i++){
        if(x[i]==isoterma){
            cout<< "encontramos la isoterma en "<< i;
        }
    }
}

// m radios y n angulos
// espero que el vector te tenga n valores
void altoHorno(int ri, int re, int m, int n, double isoterma, double ti, vector<double> te){
    // tengo que:
    // A= 1 0 0 0  segun se mjultiplique x te
    // b= (ti(x n), 0 ... 0, te)
    // x= (tm=1n=1, ....... tm=m+1n=n)
    vector<double> b= armarB(n,m,te,ti);

    // creo la matriz con los coeficientes
    vector<vector<double>> A(n*m, vector<double>(n*m));

    double c_jm1_k; 
    double c_j_k;
    double c_j1_k = (1/(((re-ri)/m)^2));
    double c_j_km1; 
    double c_j_k1;
    int fila=0;
    for(int radio=0; radio<m; radio++){
        for(int angulo=0; angulo<n; angulo++){
            vector<double>aux(n*m,0);
            if(radio==0){
                aux[angulo]=1;
                A[fila]=aux;
                fila++;
                break;
            }else if(radio==m-1){
                aux[angulo]=1;
                A[fila]=aux;
                fila++;
            }
            
            c_jm1_k= (1/(((re-ri)/m)^2) - (1/((radio*((re-ri)/m)))));
            c_j_k= (-2/(((re-ri)/m)^2))+(1/(radio *((re-ri)/m))) -(2/((2*3)/n)^2);
            c_j_km1= (1/(radio^2)*((2*3)/n)^2);
            c_j_k1= (1/(radio^2)*((2*3)/n)^2);

            // esto se hace por cada fila de la matriz

            aux[fila]=c_j_k;
            aux[fila-1]=c_j_km1;
            aux[fila+1]=c_j_k1;
            aux[fila+n]=c_j1_k;
            aux[fila-n]=c_jm1_k;

            A[fila]=aux;
            fila++;
        }
    }

    // despues deberiamos calcular la x con las factorizaciones anteriores
    vector<double> x= resolverTriangularxLU(A,b);

    encontrarIsoterma(x, isoterma);

}

int main() {
    vector<vector<double>> v = {{2,1,-1,3},{-2,0,0,0},{4,1,-2,4},{-6,-1,2,-3}};
    mostrarMatriz(facLU(v));
    vector<double> b={13,-2,24,-10};
    vector<double> x= (resolverTriangularxLU(facLU(v),b));

    //altoHorno(10,10,5,5,50,100,{20,15,18,16,25});
    return 0;
}

