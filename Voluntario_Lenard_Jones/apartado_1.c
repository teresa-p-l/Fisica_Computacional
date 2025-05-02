#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Inicializamos los parámetros

#define N 20        // Número de átomos
#define L 10.0    // Longitud de la caja
#define epsilon 1.0 
#define sigma 1.0
#define k 1.0
#define mass 1.0
#define v_0 1.0
#define h 0.002
#define Time 50
#define PI 3.14159265358979323846


// Definimos arrays

double r[N][2];             // Posiciones
double v[N][2];             // Velocidades
double a[N][2];             // Aceleraciones
double R[N][N][2];          // Distancias relativas


// Función para hacer las coordenadas periódicas

void periodic(double *r) {
    if (r[0] > L) r[0] -= L;
    if (r[0] < 0) r[0] += L;
    if (r[1] > L) r[1] -= L;
    if (r[1] < 0) r[1] += L;
}


// Establecemos las condiciones iniciales 

void initial_conditions () {
    for (int i = 0; i < N; i++) {
        r[i][0] = ((double) rand() / RAND_MAX) * L;      // Posiciones aleatorias
        r[i][1] = ((double) rand() / RAND_MAX) * L;

        double theta = ((double)rand() / (double)RAND_MAX) * 2 * PI;

        v[i][0] = v_0 * cos(theta);          // Velocidades de módulo 1 y ángulos aleatorios
        v[i][1] = v_0 * sin(theta);
    }
}


// Función para calcular la distancia relativa entre los átomos

void distance() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                R[i][j][0] = r[i][0] - r[j][0];
                R[i][j][1] = r[i][1] - r[j][1];
                periodic(R[i][j]);
                R[i][j][0] = -R[j][i][0];  
                R[i][j][1] = -R[j][i][1];
            }
        }
    }
}


// Función que calcula la aceleración por el potencial de Lennard-Jones

void aceleracion(){
    for (int i=0; i<N; i++){
        a[i][0] = 0.0;
        a[i][1] = 0.0;
        for (int j=0; j<N; j++){
            if (i != j){
                if ((R[i][j][0] > 3*sigma) || (R[i][j][1] > 3*sigma)) continue;  // Si la distancia es mayor que 3 sigma, no se considera
                else{
                    a[i][0] += 4*epsilon*(12*pow(sigma,12)/pow(R[i][j][0],13) - 6*pow(sigma,6)/pow(R[i][j][0],7))/mass;  // aquí me pone copilot un término que no le encuentro mucho sentido
                    a[i][1] += 4*epsilon*(12*pow(sigma,12)/pow(R[i][j][1],13) - 6*pow(sigma,6)/pow(R[i][j][1],7))/mass;  // aquí me pone copilot un término que no le encuentro mucho sentido
                }
            }
        }
    }
}


// Función que realiza el algoritmo de Verlet

void verlet (FILE *archivo_posiciones) {

    double omega[N][2];  // Definimos un vector auxiliar
    for (int i = 0; i < N; i++) {
        r[i][0] += v[i][0] * h + 0.5 * a[i][0] * h * h;
        r[i][1] += v[i][1] * h + 0.5 * a[i][1] * h * h;
        periodic(r[i]);
        omega[i][0] = v[i][0] + 0.5 * a[i][0] * h; 
        omega[i][1] = v[i][1] + 0.5 * a[i][1] * h;
    }

    distance();
    aceleracion();

    for (int i = 0; i < N; i++) {
        v[i][0] = omega[i][0] + 0.5 * a[i][0] * h; 
        v[i][1] = omega[i][1] + 0.5 * a[i][1] * h;
    }
}


// Hacemos una función que nos de la enerergía cinética, potencial y total y la escriba en un archivo

void compute_energy(FILE *archivo_energia) {
    double E_kin = 0.0;
    double E_pot = 0.0;
    double E_tot = 0.0;

    for (int i = 0; i < N; i++) {
        E_kin += 0.5 * mass * (v[i][0] * v[i][0] + v[i][1] * v[i][1]);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                if ((R[i][j][0] > 3*sigma) || (R[i][j][1] > 3*sigma)) continue;  // Si la distancia es mayor que 3 sigma, no se considera
                else{
                    E_pot += 4*epsilon*(pow(sigma/R[i][j][0],12) - pow(sigma/R[i][j][0],6)); // Potencial de Lennard-Jones
                }
            }
        }
    }

    E_tot = E_kin + E_pot;

    fprintf(archivo_energia, "%e %e %e\n", E_kin, E_pot, E_tot);
}




int main() {
    // tu código
    return 0;
}