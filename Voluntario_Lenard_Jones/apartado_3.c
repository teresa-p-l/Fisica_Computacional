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
#define v_0 1.0
#define h 0.002
#define Time 100
#define PI 3.14159265358979323846

#define Tmin 20     // Tiempos para el histograma de velocidades
#define Tmax 60


// Definimos arrays

double r[N][2];             // Posiciones
double v[N][2];             // Velocidades
double a[N][2];             // Aceleraciones
double Temp;             // Temperatura  
double V_promedio;  

// Función para hacer las coordenadas periódicas

void periodic(double *r) {
    if (r[0] >= L) r[0] -= L;
    if (r[0] < 0) r[0] += L;
    if (r[1] >= L) r[1] -= L;
    if (r[1] < 0) r[1] += L;
}

double suma_momentos(double *r, double *v){
    double suma_momento = 0.0;
    double sumax = 0.0;
    double sumay = 0.0;

    if (r[0] >= L ) {
        sumax += 2.0 * mass * v[0];
        //printf("%e\n", sumax);
    }
    if (r[0] <= 0 ){
        sumax += 2.0 * mass * v[0]; 
        //printf("%e\n", sumax);
    }
    if (r[1] >= L ){
        sumay += 2.0 * mass * v[1];
        //printf("%e\n", sumay);
    }
    if (r[1] <= 0 ) {
        sumay += 2.0 * mass * v[1];
        //printf("%e\n", sumay);
    } 

    return sumax + sumay;
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

void initial_conditions () {
   
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

double verlet (FILE *archivo_posiciones) {

    double omega[N][2];  // Definimos un vector auxiliar
    double momento = 0.0;
    for (int i = 0; i < N; i++) {
        r[i][0] += v[i][0] * h + 0.5 * a[i][0] * h * h;
        r[i][1] += v[i][1] * h + 0.5 * a[i][1] * h * h;
        momento += suma_momentos(r[i], v[i]);
        periodic(r[i]);
        omega[i][0] = v[i][0] + 0.5 * a[i][0] * h; 
        omega[i][1] = v[i][1] + 0.5 * a[i][1] * h;
    }

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

    return momento; // Devolvemos el momento total
}

// Función para calcular la temperatura del sistema


double compute_histogram_v_paT() {
    
    double suma = 0.0;
    for (int i = 0; i < N; i++) {
        suma += sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] );
    }
    return suma / N;
}






int main(void) {
    FILE *archivo_posiciones = fopen("posiciones.txt", "w");
    if (archivo_posiciones == NULL) {
        printf("Error al abrir el archivo de posiciones.\n");
        return 1;
    }

    double momento_final = 0.0;

    initial_conditions ();  
    aceleracion();  

 

    V_promedio = 0.0;
    // Bucle principal de la simulación
    int steps = (int)(Time/h);
    for (int step = 0; step < steps; step++) {
        momento_final += verlet(archivo_posiciones);
        if ((step >= (int)(Tmin/h)) && (step <(int)(Tmax/h))){
            V_promedio += compute_histogram_v_paT(); // Acumulamos la temperatura
        }
        
    }

    Temp = V_promedio / ((Tmax - Tmin) / h); // Temperatura final
    printf("Temperatura final: %e\n", Temp); // Imprimir la temperatura final


    printf("momento: %e\n", momento_final/L);
    
    fclose(archivo_posiciones);
    printf("Simulación completada. Resultados guardados en archivos de salida.\n");

    return 0;
}