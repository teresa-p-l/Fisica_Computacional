#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define g 9.80665
#define pi 3.14159265358979323846
#define h 0.01        // Paso temporal
#define E 10.0         // Energía total del sistema
#define tf 100        // Tiempo total de la simulacion



void condiciones_iniciales(double *y){

    y[0]=pi/16;   //Phi
    y[1]=pi/16;   //Psi
    y[2]=sqrt(-3*g+E+2*g*cos(y[0])+g*cos(y[1]));     //Velocidad de Phi 
    y[3]=0;                                       //Velocidad de Psi
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
    FILE *archivo_posiciones=fopen("posiciones_E.txt", "w");
    // Verificar que todos los archivos se abrieron correctamente
    if (archivo_posiciones == NULL) {
        printf("Error: No se pudo abrir archivo_posiciones\n");
        return 1;
    }

    FILE *archivo_phi_psi=fopen("poincare_phi_psi_E.txt", "w");
    if (archivo_phi_psi == NULL) {
        printf("Error: No se pudo abrir archivo_phi_psi\n");
        return 1;
    }

    FILE *archivo_phi_dphi=fopen("poincare_phi_dphi_E.txt", "w");
    if (archivo_phi_dphi == NULL) {
        printf("Error: No se pudo abrir archivo_phi_dphi\n");
        return 1;
    }

    FILE *archivo_psi_dpsi=fopen("poincare_psi_dpsi_E.txt", "w");
    if (archivo_psi_dpsi == NULL) {
        printf("Error: No se pudo abrir archivo_psi_dpsi\n");
        return 1;
    }




    // Definimos los vectores y variables que vamos a necesitar

    double y[4];        
    double k[4][4];              


    int i, j; //Contadores
    double t = 0; //Contador de tiempo

    
    double x1, y1;      // Posiciones de la primera masa
    double x2, y2;      // Posiciones de la segunda masa





    //Calculo e imprimo las posiciones iniciales
    condiciones_iniciales(y);

    fprintf(archivo_posiciones, "%lf, %lf \n", x1, y1);
    fprintf(archivo_posiciones, "%lf, %lf \n", x2, y2);
    fprintf(archivo_posiciones, "\n");

    //Imprimo los valores iniciales de phi y psi para hacer el mapa de Poincaré
    fprintf(archivo_phi_psi, "%lf, %lf", y[0], y[1]);
    fprintf(archivo_phi_psi, "\n");

    //Imprimo los valores iniciales de phi y dphi para hacer el mapa de Poincaré
    fprintf(archivo_phi_dphi, "%lf, %lf", y[0], y[2]);
    fprintf(archivo_phi_dphi, "\n");


 


    while(t<tf)
    {
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
            
        x1=sin(y[0]);
        y1=-cos(y[0]);
        x2=x1+sin(y[1]);
        y2=y1-cos(y[1]);

            

        //Imprimo las trayectorias en un fichero y los mapas de Poincaré en otro
        fprintf(archivo_posiciones, "%lf, %lf \n", x1, y1);
        fprintf(archivo_posiciones, "%lf, %lf \n", x2, y2);
        fprintf(archivo_posiciones, "\n");

        fprintf(archivo_phi_psi, "%lf, %lf", y[0], y[1]);
        fprintf(archivo_phi_psi, "\n");

        fprintf(archivo_phi_dphi, "%lf, %lf", y[0], y[2]);
        fprintf(archivo_phi_dphi, "\n");

        fprintf(archivo_psi_dpsi, "%lf, %lf", y[1], y[3]);
        fprintf(archivo_psi_dpsi, "\n");

            
            
      
       //Añado un paso temporal
        t+=h;
        


    }
            
    fclose(archivo_posiciones);
    fclose(archivo_phi_psi);
    fclose(archivo_phi_dphi);
    fclose(archivo_psi_dpsi);
    return 0;
}