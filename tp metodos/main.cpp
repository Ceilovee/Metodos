#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <chrono>


using namespace std;

void mostrarMatriz(vector<vector<double>> &m){
    for(int i=0;i<m.size();i++){
        for(int j=0;j<m.size();j++) {
            cout<< m[i][j]<<" ";
        }
        cout<< endl;
    }
}

vector<vector<double>> eliminacionGauss(vector<vector<double>> &m, vector<double> &b){
    for(int i=0;i<m.size();i++){
        // if es 0 mover fila
        for(int j=i+1;j<m.size();j++){
            double num =m[j][i]/m[i][i];
            for(int k=i;k<m.size()+1;k++){
                m[j][k] = m[j][k] - (num* m[i][k]);
            }
            b[j]=b[j]-(num*b[i]);
        }
        
    }
    return m;
}

vector<vector<double>> facLU(vector<vector<double>> &m){
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


vector<double> resolverTriangular_L(vector<vector<double>> &m, vector<double> &b){
    int n= m.size();
    vector<double> x(n);  
    
    for (int i = 0; i < n; i++){
        x[i] = b[i];
        
        for (int j=0; j<i; j++){
            x[i] -= m[i][j]*x[j];
        }
    }
    
    return x;
}
 
vector<double> resolverTriangular(vector<vector<double>> &m, vector<double> &b){
    int n= m.size();
    vector<double> x(n);  
    
    
    for (int i = n-1; i >= 0; i--){
        x[i] = b[i];
 
        for (int j=i+1; j<n; j++){
            x[i] -= m[i][j]*x[j];
        }
        x[i] = x[i]/m[i][i];
    }
  
    return x;
}

vector<double> resolverTriangularxLU(vector<vector<double>> &m,vector<double> &b){
    vector<double> x(b.size());
    b = resolverTriangular_L(m,b);
    x = resolverTriangular(m,b);
    return x;
}

vector<double> armarB(int m, int n, vector<double> &te){
    vector<double> b(m*n);
    for(int i=0; i<b.size(); i++){
        if(i<n){
            b[i]= te[i]; 
        }else if(i>=(b.size()-n)){
            b[i]=te[(2*n)-(b.size()-i)];
        }else{
            b[i]=0;
        }
    }
    // con esto se crea b de la forma requerida

    return b;
}

vector<vector<double>> crearA(double ri, double re, double m, double n){

    // creamos la matriz con todo en 0
    vector<vector<double>> A(n*m, vector<double>(n*m,0));
    vector<double>aux(n*m,0.0);

    // los coeficientes
    double c_jm1_k; 
    double c_j_k;
    double c_j1_k = (1/((pow(((re-ri)/(m-1)),2))));
    double c_j_km1; 
    double c_j_k1;

    int radio=0;
    for(int i=0;i<n*m;i++){
        
        radio = int( i / (n) );
       if( i < n || i>= (n*m)-(n)){
           A[i][i]=1;  
       }else {
            double radNivel= (radio*((re-ri)/(m-1)))+ri;
            c_j1_k = (1/((pow(((re-ri)/(m-1)),2))));
            c_jm1_k= (1.0/pow(((re-ri)/(m-1)),2)) - (1.0/((radNivel*((re-ri)/(m-1)))));
            c_j_k= (-2.0/((pow(((re-ri)/(m-1)),2))) + (1.0/(radNivel *((re-ri)/(m-1)))) - (2.0 / (pow(( (2.0*M_PI) / n ) , 2 )*pow(radNivel,2)) ));
            c_j_km1= (1.0/(pow(radNivel,2)*(pow((2*M_PI)/n,2))));
            c_j_k1= (1.0/(pow(radNivel,2)*(pow((2*M_PI)/n,2))));
           
           
           if ( i % int(n) == 0){

                A[i][i]=c_j_k;
                A[i][i+(n-1)]=c_j_km1;
                A[i][i+1]=c_j_k1;
                A[i][i+n]=c_j1_k;
                A[i][i-n]=c_jm1_k;

           }
           else if(i % int(n) == n-1){
                A[i][i]=c_j_k;
                A[i][i-1]=c_j_km1;
                A[i][i-(n-1)]=c_j_k1;
                A[i][i+n]=c_j1_k;
                A[i][i-n]=c_jm1_k;


           }else{
                if ( i < (n*m)-n){
                    A[i][i]=c_j_k;
                    A[i][i-1]=c_j_km1;
                    A[i][i+1]=c_j_k1;
                    A[i][i+n]=c_j1_k;
                    A[i][i-n]=c_jm1_k;
                }
           }
       }
    }
    
    return A;

}

void altoHorno(vector<vector<double>> &A, double m, double n, double isoterma, vector<double> &te, string algoritmo, string filename_out){
  
    vector<double> b= armarB( m,n,te);

    string filename(filename_out);
    fstream file_out;
    file_out.open(filename, std::ios::app);
    vector<double> x;

    if (algoritmo=="GAUSS"){
        A=eliminacionGauss(A,b);
        x= resolverTriangular(A,b);
    }
    else{
        x= resolverTriangularxLU(A,b);
    }

    // lo mandamos a un archivo
    for (int i=0; i<n*m; i++) file_out << fixed <<x[i] << endl;
    
    file_out.close();


}


int main(int argc, char** argv) {
    
    string algoritmo = argv[1];

    double ri; double re; double m; double n; double isoterma; int nist;
    std::cin >> ri >> re >> m >> n >> isoterma >> nist;
    vector<vector<double>> tes(nist,vector<double>(2*n,0));

    for(int i=0; i<nist;i++) for (int j = 0; j < 2*n; j++) cin >> tes[i][j];

    if (algoritmo=="GAUSS"){
        auto start = chrono::steady_clock::now();
        
        vector<vector<double>> A = crearA(ri,re,m,n);
        for(int i=0;i<nist;i++) altoHorno(A,m,n,isoterma,tes[i],algoritmo,argv[2]);

        auto end = chrono::steady_clock::now();
        double total_time = chrono::duration<double>(end - start).count();
        cout << total_time << endl;
    }

    else if(algoritmo=="LU"){
        auto start = chrono::steady_clock::now();
        vector<vector<double>> A = crearA(ri,re,m,n);
        A= facLU(A);
        for(int i=0;i<nist;i++) altoHorno(A,m,n,isoterma,tes[i],algoritmo, argv[2]);
        auto end = chrono::steady_clock::now();
        double total_time = chrono::duration<double>(end - start).count();
        cout << total_time << endl;
    }

    else{
        cout << "Algoritmo invalido" << endl;
        return -1;
    }
    
    return 0;
}
