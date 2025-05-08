#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Inicializamos los parámetros

#define N 16        // Número de átomos
#define L 10.0    // Longitud de la caja
#define epsilon 1.0 
#define sigma 1.0
#define k 1.0
#define mass 1.0
#define h 0.002
#define Time 100
#define PI 3.14159265358979323846
#define nmom 10


#define Tmin 20     // Tiempos para el histograma de velocidades
#define Tmax 60


// Definimos arrays

double r[N][2];             // Posiciones
double v[N][2];             // Velocidades
double a[N][2];             // Aceleraciones
double Temp;             // Temperatura  
double V_promedio;  

// Función para hacer las coordenadas periódicas
// Pasamos momento como puntero para modificar su valor

void periodic(double r[N][2], double v[N][2], double *momento) {
    for (int i = 0; i < N; i++) {
        if (r[i][0] >= L) {
            r[i][0] -= L; 
            if(v[i][0] > 0) {
                *momento += mass * v[i][0];
            } else {
                *momento -= mass * v[i][0];
            }
        }
        if (r[i][0] < 0) {
            r[i][0] += L;
            if(v[i][0] > 0) {
                *momento += mass * v[i][0];
            } else {
                *momento -= mass * v[i][0];
            }
        }
        if (r[i][1] >= L) {
            r[i][1] -= L; 
            if(v[i][1] > 0) {
                *momento += mass * v[i][1];
            } else {
                *momento -= mass * v[i][1];
            }
        }
        if (r[i][1] < 0) {
            r[i][1] += L; 
            if(v[i][1] > 0) {
                *momento += mass * v[i][1];
            } else {
                *momento -= mass * v[i][1];
            }

        }
    }
}



// Función para calcular la distancia mínima entre dos puntos con condiciones periódicas

void dist_min(double *r_i, double *r_j, double *R) {
    R[0] = r_i[0] - r_j[0];
    R[1] = r_i[1] - r_j[1];
    
    // Aplicar condiciones periódicas a la distancia relativa
    if (R[0] > L/2) R[0] -= L;
    if (R[0] < -L/2) R[0] += L;
    if (R[1] > L/2) R[1] -= L;
    if (R[1] < -L/2) R[1] += L;
}

// Establecemos las condiciones iniciales 

void initial_conditions (double v_0) {
   
    for(int i = 0; i < N; i++) {
        r[i][0] = (i % ((int)sqrt(N)+1))*L/((int)sqrt(N)+1);      // Posiciones aleatorias
        r[i][1] = (i / ((int)sqrt(N)+1))*L/((int)sqrt(N)+1); 

        double theta = ((double)rand() / (double)RAND_MAX) * 2 * PI;

        v[i][0] = v_0 * cos(theta);          // Velocidades de módulo 1 y ángulos aleatorios
        v[i][1] = v_0 * sin(theta);
    }
  
 /*
    srand(time(NULL));  // Inicializar la semilla para números aleatorios

    for (int i = 0; i < N; i++) {
        r[i][0] = ((double) rand() / RAND_MAX) * L;      // Posiciones aleatorias
        r[i][1] = ((double) rand() / RAND_MAX) * L;

        double theta = ((double)rand() / (double)RAND_MAX) * 2 * PI;

        v[i][0] = v_0 * cos(theta);          // Velocidades de módulo 1 y ángulos aleatorios
        v[i][1] = v_0 * sin(theta);
    }
*/
        
}


// Función que calcula la aceleración por el potencial de Lennard-Jones

void aceleracion(){
    for (int i=0; i<N; i++){    // Inicializamos la aceleración a cero
        a[i][0] = 0.0;
        a[i][1] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {  // Solo calcular cada par una vez

            double R[2];    // Distancia relativa
            dist_min( r[i], r[j], R );

            double mod_R = sqrt( R[0]*R[0] + R[1]*R[1] );
            
            
            if ( R[0]*R[0] + R[1]*R[1] < 9.0*sigma*sigma ) { // Solo considerar partículas dentro del rango de interacción (r < 3*sigma)

                // Calculamos la fuerza sin dirección
                double f_mag = 4.0 * epsilon * ( 12.0 * pow(sigma/mod_R, 12) - 6.0 * pow(sigma/mod_R, 6)) / mod_R;
                
                // Componentes de la fuerza
                double fx = f_mag * R[0] / mod_R;
                double fy = f_mag * R[1] / mod_R;
                
                // Aplicamos la segunda ley de Newton
                a[i][0] += fx / mass;
                a[i][1] += fy / mass;
                a[j][0] -= fx / mass;
                a[j][1] -= fy / mass;
            }
        }
    }
}


// Función que realiza el algoritmo de Verlet
// Ahora recibe un puntero a momento

void verlet(FILE *archivo_posiciones, double *momento) {
    double omega[N][2];  // Definimos un vector auxiliar
    for (int i = 0; i < N; i++) {
        r[i][0] += v[i][0] * h + 0.5 * a[i][0] * h * h;
        r[i][1] += v[i][1] * h + 0.5 * a[i][1] * h * h;
        
        omega[i][0] = v[i][0] + 0.5 * a[i][0] * h; 
        omega[i][1] = v[i][1] + 0.5 * a[i][1] * h;
    }
    
    periodic(r, v, momento);
    aceleracion();

    for (int i = 0; i < N; i++) {
        v[i][0] = omega[i][0] + 0.5 * a[i][0] * h; 
        v[i][1] = omega[i][1] + 0.5 * a[i][1] * h;
    }

    // Añadimos las posiciones al archivo de salida
    for (int i = 0; i < N; i++) {
        fprintf(archivo_posiciones, "%e %e\n", r[i][0], r[i][1]);
    }
    fprintf(archivo_posiciones, "\n");
}

// Función para calcular la temperatura del sistema

double compute_histogram_v_paT() {
    double suma = 0.0;
    for (int i = 0; i < N; i++) {
        suma += sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1]);
    }
    return suma / N;
}




void compute_energy(FILE *archivo_energia) {
    double E_kin = 0.0;
    double E_pot = 0.0;
    double E_tot=0.0;

    for (int i = 0; i < N; i++) {
        E_kin += 0.5 * mass * (v[i][0] * v[i][0] + v[i][1] * v[i][1]);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {  // Solo calcular cada par una vez
            
            if (i != j) {
                double R[2];
                dist_min(r[i], r[j], R);
            
                double mod_R = sqrt( R[0]*R[0] + R[1]*R[1] );
            
                if ( R[0]*R[0] + R[1]*R[1] < 9.0*sigma*sigma ) {  // Solo para r < 3*sigma
                    E_pot += 4.0 * epsilon * (pow(sigma/mod_R, 12) - pow(sigma/mod_R, 6));
                }
            
            }
        }
    }
    
    E_tot = E_kin + E_pot;

    fprintf(archivo_energia, "%e %e %e\n", E_kin, E_pot, E_tot);
}


int main(void) {
    FILE *archivo_posiciones = fopen("posiciones.txt", "w");
    if (archivo_posiciones == NULL) {
        printf("Error al abrir el archivo de posiciones.\n");
        return 1;
    }
    FILE *archivo_momentos = fopen("momentos.txt", "w");
    if (archivo_momentos == NULL) {
        printf("Error al abrir el archivo de momentos.\n");
        return 1;
    }
    FILE *archivo_energia = fopen("energia.txt", "w");
    if (archivo_energia == NULL) {
        printf("Error al abrir el archivo de energía.\n");
        return 1;
    }
    
    double v_0 = 0.0;
    fprintf(archivo_momentos, "Velocidad Temperatura Presion \n");

    for (int m = 0; m < nmom; m++) {
        v_0++;

        initial_conditions(v_0);  
        aceleracion();  
        
        double momento_total = 0.0;  // Inicializamos el momento para cada v_0
        V_promedio = 0.0;
        
        // Bucle principal de la simulación
        int steps = (int)((Time-Tmin)/h);
        for (int step = 0; step < steps; step++) {
            verlet(archivo_posiciones, &momento_total);  // Pasamos el puntero a momento_total
            compute_energy(archivo_energia);
            V_promedio += compute_histogram_v_paT(); // Acumulamos la temperatura
        }
    
        Temp = V_promedio*V_promedio/2 / ((Time-Tmin) / h); // Temperatura final
        fprintf(archivo_momentos, "%f %f %e \n", v_0, Temp, momento_total/(L*((Time-Tmin)/h))); // Imprimir la temperatura final
    }

    FILE *datos_simulacion = fopen("datos_simulacion.txt", "w");
    fprintf(datos_simulacion, "%d %f %f %f \n", N, L, h, Temp); // Guardar los parámetros de la simulación
    fclose(datos_simulacion);
    
    fclose(archivo_posiciones);
    fclose(archivo_momentos);
    printf("Simulación completada. Resultados guardados en archivos de salida.\n");

    return 0;
}