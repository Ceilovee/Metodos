#include <iostream>
#include <fstream>
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

vector<vector<double>> eliminacionGauss(vector<vector<double>> m, vector<double> &b){
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
  
    return x;
}

vector<double> resolverTriangularxLU(vector<vector<double>> m,vector<double> b){
    vector<double> x(b.size());
    b = resolverTriangular_L(m,b);
    x = resolverTriangular(m,b);
       
    //printf("\nSolution for the system (L):\n");
    //for (int i=0; i<x.size(); i++)
     //   printf("%lf\n", x[i]);
    return x;
}

vector<double> armarB(int n, int m, vector<double> te){
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

void encontrarIsoterma(vector<double> x, double isoterma){
    // hacer que busque bien la isoterma
    for(int i=0; i < x.size();i++){
        if(x[i]==isoterma){
            cout<< "encontramos la isoterma en "<< i;
        }
    }
}

vector<vector<double>> crearA(double ri, double re, double m, double n){

    // creo la matriz con los coeficientes
    vector<vector<double>> A(n*m, vector<double>(n*m));
    vector<double>aux(n*m,0.0);

    // los coeficientes
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
            c_j_k= (-2.0/((pow(((re-ri)/m),2))) + (1/(radNivel *((re-ri)/m))) - (2.0 / (pow(( (2.0*M_PI) / n ) , 2 )*pow(radNivel,2)) ));
            c_j_km1= (1.0/(pow(radNivel,2)*(pow((2*M_PI)/n,2))));
            c_j_k1= (1.0/(pow(radNivel,2)*(pow((2*M_PI)/n,2))));

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
    return A;

}

// m radios y n angulos
// espero que el vector te tenga n valores
void altoHorno(vector<vector<double>> A, int m, int n, double isoterma, vector<double> te){
    // tengo que:
    // A= 1 0 0 0  segun se mjultiplique x te
    // b= (ti(x n), 0 ... 0, te)
    // x= (tm=1n=1, ....... tm=m+1n=n)
    //cout << "hola" <<endl;
    vector<double> b= armarB(n,m,te);

    string filename("testResultado.txt");
    fstream file_out;
    file_out.open(filename, std::ios::out);// | std::ios::app);
    A=eliminacionGauss(A,b);
    vector<double> x= resolverTriangular(A,b);
    // lo mandamos a un archivo
    for (int i=0; i<n*m; i++) file_out << fixed <<x[i] << endl;
    
    file_out.close();

    //encontrarIsoterma(x, isoterma);

}

int main() {
    
    vector<vector<double>> v = {{2,1,-1,3},{-2,0,0,0},{4,1,-2,4},{-6,-1,2,-3}};
    mostrarMatriz(facLU(v));
    vector<double> b={13,-2,24,-10};
    v= eliminacionGauss(v,b);
    for(int i=0;i<b.size();i++) cout<< b[i]<<",";
    vector<double> x= resolverTriangular(v,b);
    for(int i=0;i<x.size();i++) cout<< x[i]<<",";
    
    //solucion {1,-30,7,16}
    
    //altoHorno(10,20,5,5,50,100,{20,15,18,16,25});
    double ri; double re; double m; double n; double isoterma; int nist;
    std::cin >> ri >> re >> n >> m >> isoterma >> nist;
    vector<vector<double>> tes(nist,vector<double>(2*n,0));
    // aca cambie a 2*n porque tiene a te y ti
    for(int i=0; i<nist;i++) for (int j = 0; j < 2*n; j++) cin >> tes[i][j];


    //vector<vector<double>> A = crearA(10,100,30,30);
    vector<vector<double>> A = crearA(ri,re,m,n);
    //A= facLU(A);
    //nist=1;
    //for(int i=0;i<nist;i++) altoHorno(A,30,30,isoterma,{1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0});
    for(int i=0;i<nist;i++) altoHorno(A,m,n,isoterma,tes[i]);
    return 0;
}
