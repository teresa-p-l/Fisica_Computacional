#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define g 9.80665
#define pi 3.14159265358979323846
#define h 0.01        // Paso temporal
#define tf 12000        // Tiempo total de la simulacion
#define perturbacion 0.000000000000001 // Perturbación inicial




void condiciones_iniciales(double *y, double *ypert, double E){

        y[0]=0.1;   //Phi
        y[1]=0.6;   //Psi
        y[2]=sqrt(-3*g+E+2*g*cos(y[0])+g*cos(y[1]));     //Velocidad de Phi 
        y[3]=0;                                     //Velocidad de Psi

        ypert[0]=0.1;   //Psi
        ypert[1]=0.601;   //Psi         //Psi
        ypert[2]=sqrt(-3*g+E+2*g*cos(y[0])+g*cos(y[1])) + 0.00000;         //Velocidad de Phi
        ypert[3]= 0.00000;         //Velocidad de Psi
}

// Hacemos 4 funciones para calcular las derivadas primera y segunda de phi y psi para que el Runge-Kutta sea más sencillo

double f_dphi(double dphi)   // Calcula la primera derivada de phi
{
    return dphi;
}

double f_dpsi(double dpsi)  // Calcula la primera derivada de psi
{
    return dpsi;
}


double f_ddphi(double phi, double psi, double dphi, double dpsi)    // Calcula la segunda derivada de phi
{
    double ddphi= (g*sin(psi)*cos(phi-psi) - 2*g*sin(phi) - pow(dphi, 2)*cos(phi-psi)*sin(phi-psi) - pow(dpsi, 2)*sin(phi - psi))/(2 - pow(cos(phi-psi), 2));
    
    return ddphi;
}

double f_ddpsi(double phi, double psi, double dphi, double dpsi)    // Calcula la segunda derivada de psi
{
    double ddpsi = (g*sin(phi)*cos(phi-psi) - g*sin(psi) + 0.5*pow(dpsi, 2)*cos(phi-psi)*sin(phi-psi) + pow(dphi, 2)*sin(phi-psi))/(1 - 0.5*pow(cos(phi-psi), 2));
  
    return ddpsi;
}







int main(void)
{

double E=15.0;




    char nombre_archivo_distancia[50]; // Buffer para el nombre del archivo
    sprintf(nombre_archivo_distancia, "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/distancia_%.0f.txt",E); // Construir el nombre dinámicamente

    FILE *archivo_distancia=fopen(nombre_archivo_distancia, "w");
    if (archivo_distancia == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }

    char nombre_archivo_lambda[50]; // Buffer para el nombre del archivo
    sprintf(nombre_archivo_lambda, "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_%.0f.txt",E); // Construir el nombre dinámicamente

    FILE *archivo_lambda=fopen(nombre_archivo_lambda, "w");
    if (archivo_lambda == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }

    char nombre_archivo_t[50]; // Buffer para el nombre del archivo
    sprintf(nombre_archivo_t, "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/t_%.0f.txt",E); // Construir el nombre dinámicamente

    FILE *archivo_t=fopen(nombre_archivo_t, "w");
    if (archivo_t == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }



    // Definimos los vectores y variables que vamos a necesitar

    double y[4];        
    double k[4][4];
            
    double ypert[4];
    double kpert[4][4];
    double lambda = 0;

    int i, j; //Contadores
    double t = 0; //Contador de tiempo

            
    double x1, y1;      // Posiciones de la primera masa
    double x2, y2;      // Posiciones de la segunda masa
    double lyapu = 0.0;




    //Calculo e imprimo las posiciones iniciales
    condiciones_iniciales(y, ypert, E);

    double delta = sqrt( pow((y[0]-ypert[0]),2) + pow((y[1]-ypert[1]),2) + pow((y[2]-ypert[2]),2) + pow((y[3]-ypert[3]),2) );
    double delta0 = sqrt( pow((y[0]-ypert[0]),2) + pow((y[1]-ypert[1]),2) + pow((y[2]-ypert[2]),2) + pow((y[3]-ypert[3]),2) );

    fprintf(archivo_distancia, "%lf \n", delta);


    for (int tiempototal = 100; tiempototal < 1000; tiempototal+=100) {

        t=0;
        delta = 0.0;
        lyapu = 0.0;
       while(t<tiempototal) {    //  CAMBIAR AQUÍ TF POR TIEMPOTOTAL
        

            //Calculo de k1:
            k[0][0]=h*f_dphi(y[2]);
            k[0][1]=h*f_dpsi(y[3]);
            k[0][2]=h*f_ddphi(y[0], y[1], y[2], y[3]);
            k[0][3]=h*f_ddpsi(y[0], y[1], y[2], y[3]);
        
            //Calculo de k2:
            k[1][0]=h*f_dphi(y[2]+0.5*k[0][2]);
            k[1][1]=h*f_dpsi(y[3]+0.5*k[0][3]);
            k[1][2]=h*f_ddphi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            k[1][3]=h*f_ddpsi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
        
            //Calculo de k3:
            k[2][0]=h*f_dphi(y[2]+0.5*k[1][2]);
            k[2][1]=h*f_dpsi(y[3]+0.5*k[1][3]);
            k[2][2]=h*f_ddphi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);
            k[2][3]=h*f_ddpsi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);
            //Calculo de k4
            k[3][0]=h*f_dphi(y[2]+k[2][2]);
            k[3][1]=h*f_dpsi(y[3]+k[2][3]);
            k[3][2]=h*f_ddphi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            k[3][3]=h*f_ddpsi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);

            //Calculo de la nueva y
            for(i=0;i<4;i++) {
                y[i]=y[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;
            }
            




            // PARA EL PERTURBADO
            //Calculo de k1:
            kpert[0][0]=h*f_dphi(ypert[2]);
            kpert[0][1]=h*f_dpsi(ypert[3]);
            kpert[0][2]=h*f_ddphi(ypert[0], ypert[1], ypert[2], ypert[3]);
            kpert[0][3]=h*f_ddpsi(ypert[0], ypert[1], ypert[2], ypert[3]);
            
            //Calculo de kpert2:
            kpert[1][0]=h*f_dphi(ypert[2]+0.5*kpert[0][2]);
            kpert[1][1]=h*f_dpsi(ypert[3]+0.5*kpert[0][3]);
            kpert[1][2]=h*f_ddphi(ypert[0]+0.5*kpert[0][0], ypert[1]+0.5*kpert[0][1], ypert[2]+0.5*kpert[0][2], ypert[3]+0.5*kpert[0][3]);
            kpert[1][3]=h*f_ddpsi(ypert[0]+0.5*kpert[0][0], ypert[1]+0.5*kpert[0][1], ypert[2]+0.5*kpert[0][2], ypert[3]+0.5*kpert[0][3]);
        
            //Calculo de kpert3:
            kpert[2][0]=h*f_dphi(ypert[2]+0.5*kpert[1][2]);
            kpert[2][1]=h*f_dpsi(ypert[3]+0.5*kpert[1][3]);
            kpert[2][2]=h*f_ddphi(ypert[0]+0.5*kpert[1][0], ypert[1]+0.5*kpert[1][1], ypert[2]+0.5*kpert[1][2], ypert[3]+0.5*kpert[1][3]);
            kpert[2][3]=h*f_ddpsi(ypert[0]+0.5*kpert[1][0], ypert[1]+0.5*kpert[1][1], ypert[2]+0.5*kpert[1][2], ypert[3]+0.5*kpert[1][3]);
            //Calculo de kpert4
            kpert[3][0]=h*f_dphi(ypert[2]+kpert[2][2]);
            kpert[3][1]=h*f_dpsi(ypert[3]+kpert[2][3]);
            kpert[3][2]=h*f_ddphi(ypert[0]+kpert[2][0], ypert[1]+kpert[2][1], ypert[2]+kpert[2][2], ypert[3]+kpert[2][3]);
            kpert[3][3]=h*f_ddpsi(ypert[0]+kpert[2][0], ypert[1]+kpert[2][1], ypert[2]+kpert[2][2], ypert[3]+kpert[2][3]);

            //Calculo de la nueva ypert
            for(i=0;i<4;i++) {
                ypert[i]=ypert[i]+(kpert[0][i]+2*kpert[1][i]+2*kpert[2][i]+kpert[3][i])/6;
            }




                // Calulamos la distancia entre las trayectorias y las guardamos
            delta += sqrt( pow((y[0]-ypert[0]),2) + pow((y[1]-ypert[1]),2) + pow((y[2]-ypert[2]),2) + pow((y[3]-ypert[3]),2) );
            //delta += sqrt( pow((y[0]-ypert[0]),2)  + pow((y[2]-ypert[2]),2) );
            lyapu += log(delta/delta0);
            //fprintf(archivo_lambda, "%lf \n", log(delta/delta0)/t*1.0);

            
            t+=h;
            
            



        }

    
    fprintf(archivo_distancia, "%lf \n", delta/tiempototal*1.0);
    fprintf(archivo_lambda, "%lf \n", lyapu/tiempototal*1.0);
    }

        
    fclose(archivo_distancia);

    return 0;
}