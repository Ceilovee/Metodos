#include <iostream>
#include <vector>
#include <cmath>
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
 
    printf("\nSolution for the system (L):\n");
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
    //mostrarMatriz(m);
 
    printf("\nSolution for the system:\n");
    for (int i=0; i<n; i++)
        printf("%lf\n", x[i]);
    
    return x;
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
    // hacer que busque bien la isoterma
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
    vector<double>aux(n*m,0.0);

    double c_jm1_k; 
    double c_j_k;
    double c_j1_k = (1/((pow(((re-ri)/m),2))));
    double c_j_km1; 
    double c_j_k1;
    int fila=0;
    for(int radio=0; radio<m; radio++){
        for(int angulo=0; angulo<n; angulo++){
            for(int k=0; k< aux.size();k++){
                aux[k]=0.0;
            }
            if(radio==0){
                aux[angulo]=1;
                A[fila]=aux;
                fila++;
            }else if(radio==m-1){
                aux[angulo+(radio*n)]=1;
                A[fila]=aux;
                fila++;
            }
            else{
            double radNivel= (radio*((re-ri)/m))+ri;
            c_jm1_k= (1.0/pow(((re-ri)/m),2)) - (1/((radNivel*((re-ri)/m))));
            c_j_k= (-2.0/((pow(((re-ri)/m),2))) + (1/(radNivel *((re-ri)/m))) - ( 2.0 / (pow(( (2.0*M_PI) / n ) , 2 )*pow(radNivel,2))));
            c_j_km1= (1.0/(pow(radNivel,2)*(pow((2*M_PI)/n,2))));
            c_j_k1= (1.0/(pow(radNivel,2) * pow( ((2*M_PI)/n) ,2)));

            // esto se hace por cada fila de la matriz

            aux[fila]=c_j_k;
            aux[fila-1]=c_j_km1;
            aux[fila+1]=c_j_k1;
            aux[fila+n]=c_j1_k;
            aux[fila-n]=c_jm1_k;

            A[fila]=aux;
            fila=fila+1;
            }
            
        }
    }

    // despues deberiamos calcular la x con las factorizaciones anteriores
     vector<double> x= resolverTriangularxLU(facLU(A),b);

     encontrarIsoterma(x, isoterma);

}

int main() {
    /*vector<vector<double>> v = {{2,1,-1,3},{-2,0,0,0},{4,1,-2,4},{-6,-1,2,-3}};
    mostrarMatriz(facLU(v));
    vector<double> b={13,-2,24,-10};
    vector<double> x= (resolverTriangularxLU(facLU(v),b));
    //solucion {1,-30,7,16}*/
    
    //altoHorno(10,20,5,5,50,100,{20,15,18,16,25});
    int ri; int re; int m; int n; double isoterma; double ti;
    std::cin >> ri >> re >> n >> m >> isoterma >> ti;
    vector<double> te(n,0);
    for (int i = 0; i < n; ++i) cin >> te[i];

    altoHorno(ri,re,m,n,isoterma,ti,te);
    return 0;
}
