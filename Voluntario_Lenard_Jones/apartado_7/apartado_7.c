#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Inicializamos los parámetros

#define N 16        // Número de átomos
#define L 4.0    // Longitud de la caja
#define epsilon 1.0 
#define sigma 1.0
#define k 1.0
#define mass 1.0
#define h 0.002
#define Time 250.0
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
   
    // Inicializamos las posiciones y velocidades de las partículas
    // Se pueden cambiar las condiciones iniciales aquí para probar diferentes configuraciones
    

    for(int i = 0; i < N; i++) {
        r[i][0] = (i % ((int)sqrt(N))) * ( L / ((int)sqrt(N))) + 0.5;      // Posiciones aleatorias
        r[i][1] = (i / ((int)sqrt(N))) * ( L / ((int)sqrt(N))) + 0.5; 

        double theta = ((double)rand() / (double)RAND_MAX) * 2 * PI;

        v[i][0] = v_0 * cos(theta);          // Velocidades de módulo 1 y ángulos aleatorios
        v[i][1] = v_0 * sin(theta);
    }

}


// Función que calcula la aceleración por el potencial de Lennard-Jones

void aceleracion(){
    for (int i=0; i<N; i++){    // Inicializamos la aceleración a cero
        a[i][0] = 0.0;
        a[i][1] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {  // Solo calcular cada par una vez
            if (i != j){
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


// Hacemos una función que nos de la enerergía cinética, potencial y total y la escriba en un archivo

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


// Función para calcular la temperatura del sistema y el histograma de velocidades


double compute_histogram_v_paT() {
    
    double suma = 0.0;
    for (int i = 0; i < N; i++) {
        suma += sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] );
    }
    return suma / N;
}


void compute_histogram_v(FILE *archivo_velocidades) {
    for (int i = 0; i < N; i++) {
        fprintf(archivo_velocidades, "%e\n", sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] ));
    }

}


void reescalar_v(double factor){
    for (int i = 0; i < N; i++) {
        v[i][0] *= factor;
        v[i][1] *= factor;
    }
}

// Función para calcular el cuadrado del desplazamiento respecto a su posición inicial

double displacement(FILE *archivo_desplazamiento, int i){
    double dx, dy, d = 0.0;

                dx = r[i][0] - r[0][0];
                dy = r[i][1] - r[0][1];
                d = sqrt(dx*dx + dy*dy);    
    return(d); // Guardamos el desplazamiento medio al cuadrado
}

double probemos(){
    double probar = 0.0;
    for (int i = 0; i < N; i++) {
        probar += 0.5 * mass * (v[i][0] * v[i][0] + v[i][1] * v[i][1]);
    }

    return probar/(N*1.0);
}












int main(void) {
    FILE *archivo_posiciones = fopen("posiciones.txt", "w");
    if (archivo_posiciones == NULL) {
        printf("Error al abrir el archivo de posiciones.\n");
        return 1;
    }

    FILE *archivo_energia = fopen("energia.txt", "w");
    if (archivo_energia == NULL) {
        printf("Error al abrir el archivo de energía.\n");
        return 1;
    }

    FILE *archivo_velocidades = fopen("velocidades.txt", "w");
    if (archivo_velocidades == NULL) {
        printf("Error al abrir el archivo de velocidades.\n");
        return 1;
    }

    FILE *archivo_desplazamiento = fopen("desplazamiento.txt", "w");
    if (archivo_desplazamiento == NULL) {
        printf("Error al abrir el archivo de desplazamiento.\n");
        return 1;
    }



    double v_0 = 0.0;
    initial_conditions (v_0); 
    double paprobar = 0.0; 
    double r_0[N][2]; // Guardamos las posiciones iniciales

    for (int i = 0; i < N; i++) {
        r_0[i][0] = r[i][0];
        r_0[i][1] = r[i][1]; 
        fprintf(archivo_posiciones, "%e %e\n", r[i][0], r[i][1]);
    }
    fprintf(archivo_posiciones, "\n");

    aceleracion();  

    

    V_promedio = 0.0;
    double desplazamiento_medio = 0.0;
    // Bucle principal de la simulación
    int steps = (int)(Time/h);
    for (int step = 0; step < steps; step++) {
        verlet(archivo_posiciones);
        compute_energy(archivo_energia);
        if ((step == 20/h) || (step == 40/h) || (step == 60/h) || (step == 80/h) || (step == 100/h) || (step == 120/h) || (step == 140/h) || (step == 160/h) || (step == 180/h) || (step == 200/h) || (step == 220/h) || (step == 240/h)) reescalar_v(1.1);
   
       // if ((step >= (int)(Tmin/h)) && (step <(int)(Tmax/h))){
            compute_histogram_v(archivo_velocidades); 
            V_promedio += compute_histogram_v_paT(); // Acumulamos la temperatura
            paprobar += probemos(); // Acumulamos la temperatura
        //}


        for (int i = 0; i < N; i++) {
            desplazamiento_medio += displacement(archivo_desplazamiento, i); // Acumulamos el desplazamiento medio al cuadrado
            if (step % 5 ==0){
                fprintf(archivo_desplazamiento, "%e\n", desplazamiento_medio/5.0); 
                desplazamiento_medio = 0.0; // Reiniciamos el desplazamiento medio
            }
        }
        if (step % 5 ==0) fprintf(archivo_desplazamiento, "\n");
        



        

    }


 

    Temp = paprobar / (double)(steps); 
    //Temp =  V_promedio / (double)(steps); // Temperatura final

    printf("Temperatura final: %e\n", Temp ); // Imprimir la temperatura final


    FILE *datos_simulacion = fopen("datos_simulacion.txt", "w");
    fprintf(datos_simulacion, "%d %f %f %f \n",N, L, h, Temp); // Guardar los parámetros de la simulación
    fclose(datos_simulacion);

    
    fclose(archivo_posiciones);
    fclose(archivo_energia);
    fclose(archivo_velocidades);
    fclose(archivo_desplazamiento);
    printf("Simulación completada. Resultados guardados en archivos de salida.\n");

    return 0;
}