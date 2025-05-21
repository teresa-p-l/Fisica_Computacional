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
#define Time 60.0
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



void inicializar_posiciones_hexagonal(double r[][2]) {
    // Factor de escala para ajustar a la caja de tamaño L
    double escala = L / 4.0;
    
    // Parámetros para crear la red hexagonal
    double hex_size = escala;  // Tamaño del hexágono
    
    // Constantes para la geometría hexagonal
    double dx = hex_size;      // Distancia horizontal entre centros
    double dy = hex_size * sqrt(3); // Distancia vertical entre centros
    
    // Posiciones de las partículas que forman los vértices de los hexágonos
    // Esta estructura crea un patrón de panal de abeja completo
    double posiciones[12][2] = {
        // Primer hexágono
        {1.0, 0.0},            // 0
        {0.5, 0.866},          // 1 
        {1.5, 0.866},          // 2
        
        // Segundo hexágono (compartiendo vértices con el primero)
        {2.0, 0.0},            // 3
        {2.5, 0.866},          // 4
        
        // Tercer hexágono (debajo del primero)
        {0.5, 2.598},          // 5
        {1.5, 2.598},          // 6
        
        // Cuarto hexágono (compartiendo vértices con el tercero)
        {2.5, 2.598},          // 7
        
        // Completando la estructura periódica
        {0.0, 1.732},          // 8
        {3.0, 1.732},          // 9
        {1.0, 3.464},          // 10
        {2.0, 3.464}           // 11
    };

    // Escalamos y ajustamos las posiciones a la caja
    for (int i = 0; i < N; i++) {
        // Escalamos las coordenadas
        r[i][0] = posiciones[i][0] * escala;
        r[i][1] = posiciones[i][1] * escala;
        
        // Aplicamos condiciones periódicas
        r[i][0] = fmod(r[i][0], L);
        r[i][1] = fmod(r[i][1], L);
        
        // Ajustamos valores negativos
        if (r[i][0] < 0) r[i][0] += L;
        if (r[i][1] < 0) r[i][1] += L;
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


double probemos(){
    double probar = 0.0;
    for (int i = 0; i < N; i++) {
        probar += 0.5 * mass * (v[i][0] * v[i][0] + v[i][1] * v[i][1]);
    }

    return probar/(N*1.0);
}












int main(void) {
    FILE *archivo_posiciones = fopen("posicioneshex.txt", "w");
    if (archivo_posiciones == NULL) {
        printf("Error al abrir el archivo de posiciones.\n");
        return 1;
    }

    FILE *archivo_energia = fopen("energiahex.txt", "w");
    if (archivo_energia == NULL) {
        printf("Error al abrir el archivo de energía.\n");
        return 1;
    }

    FILE *archivo_velocidades = fopen("velocidadeshex.txt", "w");
    if (archivo_velocidades == NULL) {
        printf("Error al abrir el archivo de velocidades.\n");
        return 1;
    }




    double v_0 = 0.0;
    inicializar_posiciones_hexagonal (r); 
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
    // Bucle principal de la simulación
    int steps = (int)(Time/h);
    for (int step = 0; step < steps; step++) {
        verlet(archivo_posiciones);
        compute_energy(archivo_energia);

       //if ((step >= (int)(Tmin/h)) && (step <(int)(Tmax/h))){
            compute_histogram_v(archivo_velocidades); 
            V_promedio += compute_histogram_v_paT(); // Acumulamos la temperatura
            paprobar += probemos(); // Acumulamos la temperatura
        //}
    }


 

    Temp = paprobar / (double)(steps); 
    //Temp =  V_promedio / (double)(steps); // Temperatura final

    printf("Temperatura final: %e\n", Temp ); // Imprimir la temperatura final


    FILE *datos_simulacion = fopen("datos_simulacionhex.txt", "w");
    fprintf(datos_simulacion, "%d %f %f %f \n",N, L, h, Temp); // Guardar los parámetros de la simulación
    fclose(datos_simulacion);

    
    fclose(archivo_posiciones);
    fclose(archivo_energia);
    fclose(archivo_velocidades);
    printf("Simulación completada. Resultados guardados en archivos de salida.\n");

    return 0;
}