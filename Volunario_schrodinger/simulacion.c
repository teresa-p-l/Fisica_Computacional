#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>


// Inicializamos los parámetros

#define N 110        
#define nc N/4 //número de ciclos  
#define h 1/N 
#define lamba 1.0
#define n 20
#define T 250

#define pi 3.14159265358979323846





int main (){

    double k = 2 * pi * nc / N;
    double s = 1 / ( 4*k*k);
    double V[N];
    complex phi[N][n];
    complex alpha[N-1];
    complex beta[N-1];

    for (int j=0; j<N; j++){
        if (2*N*1.0/5.0 < j && j < 3*N*1.0/5.0){
            V[j]= lamba*k*k;  
        }
        else{
            V[j]= 0;
        }
        phi[j][0] = cexp(I * k * j) * exp(-8* pow((4*j-N),2) / (N*N));
    }
    phi[0][0] = 0;
    phi[N-1][0] = 0;













}