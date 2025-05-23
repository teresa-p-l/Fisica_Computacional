#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define g 9.80665
#define pi 3.14159265358979323846
#define h 0.01        //Paso temporal
#define E 1.0         //Energía total del sistema



void condiciones_iniciales(double *y){

        y[0]=pi;  //Phi
        y[1]=pi;   //Psi
        y[2]=sqrt(E+2*g*cos(y[0])+g*cos(y[1]));     //Velocidad de Phi 
        y[3]=0;                                     //Velocidad de Psi
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

    FILE *archivo_posiciones=fopen("pendulodoble.txt", "w"); 
    FILE *archivo_phi_psi=fopen("poincare_phi_psi_E=15.txt", "w");
    FILE *archivo_phi_dphi=fopen("poincare_phi_phipunto_E=15.txt", "w");
    FILE *f4=fopen("lyapunov_E=15.txt", "w");

    
    //clock_t begin= clock(); //Empiezo a contar el tiempo de simulación

    // Definimos los vectores y variables que vamos a necesitar

    double y[4];        
    double k[4][4];              


    int i, j; //Contadores
    double t; //Contador de tiempo
    double tf; //Tiempo final de la simulación

    
    double x1, y1;      // Posiciones de la primera masa
    double x2, y2;      // Posiciones de la segunda masa




    //Variables para calculas los coeficientes de Lyapunov

    double ypert[4]; //Vector y perturbado
    double perturbacion[4] = {-0.05, -0.05, -0.05, -0.05}; //Vector que introduce las perturbaciones
    double sumadivergente = 0.0;
    double norma_perturbacion; //Lo utilizaremos para normalizar
    double exponente_lyapunov;
    int n; //Contador de iteraciones al calcular el exponente de Lyapunov

    

    //Establecemos las variables que regulan el tiempo de simulación
    tf=100;
    t=0;


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

  






    //Inicializo el vector y perturbado a partir de la perturbación introducida
    //for (i = 0; i < 4; i++) 
      //  ypert[i] = y[i] + perturbacion[i];
    

    //Para calcular los coeficientes de Lyapunov para diferentes tiempos dejamos este bucle temporal en el que cada vez que se recorre se inicializan las variables
    //necesarias para ccalcular los coeficientes de Lyapunov. 
    //Por el contrario, si no queremos calcularlos y ver solo la animación, dejamos solo la aplicación del algoritmo de Runge-Kutta (while) para el vector y, terminando
    //por imprimir las posiciones del péndulo e incrementando un paso temporal la variable t. (t+=h)

    for ( tf = 1000; tf <= 10000; tf+=1000)
    {
        t=0;
        sumadivergente=0;
        n=0;
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
            for(i=0;i<4;i++)
                y[i]=y[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;
            
            
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

            

            
            //Procedo a computar las trayectorias perturbadas

            //Calculo de k1:
            k[0][0]=h*f_dphi(ypert[2]);
            k[0][1]=h*f_dpsi(ypert[3]);
            k[0][2]=h*f_ddphi(ypert[0], ypert[1], ypert[2], ypert[3]);
            k[0][3]=h*f_ddpsi(ypert[0], ypert[1], ypert[2], ypert[3]);
        
            //Calculo de k2:
            k[1][0]=h*f_dphi(ypert[2]+0.5*k[0][2]);
            k[1][1]=h*f_dpsi(ypert[3]+0.5*k[0][3]);
            k[1][2]=h*f_ddphi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            k[1][3]=h*f_ddpsi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            
            //Calculo de k3:
            k[2][0]=h*f_dphi(ypert[2]+0.5*k[1][2]);
            k[2][1]=h*f_dpsi(ypert[3]+0.5*k[1][3]);
            k[2][2]=h*f_ddphi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);
            k[2][3]=h*f_ddpsi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);

            //Calculo de k4
            k[3][0]=h*f_dphi(ypert[2]+k[2][2]);
            k[3][1]=h*f_dpsi(ypert[3]+k[2][3]);
            k[3][2]=h*f_ddphi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            k[3][3]=h*f_ddpsi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            
            //Calculo de la nueva ypert
            for(i=0;i<4;i++)
                ypert[i]=ypert[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;



            for (i = 0; i < 4; i++) 
            {
                perturbacion[i] = ypert[i] - y[i]; //Diferencia entre las trayectorias
            }

            norma_perturbacion = sqrt(perturbacion[0] * perturbacion[0] + perturbacion[1] * perturbacion[1] +
                                    perturbacion[2] * perturbacion[2] + perturbacion[3] * perturbacion[3]);
            sumadivergente += log(fabs(norma_perturbacion)/ 0.05);

            for (i = 0; i < 4; i++) 
            {
                perturbacion[i] *= 0.05/ norma_perturbacion;
                ypert[i] = y[i] + perturbacion[i];
            }

            n++;
      
           //Añado un paso temporal
            t+=h;
        }

        exponente_lyapunov=sumadivergente/(n*h) ;
        fprintf(f4,"%lf, %lf",tf, exponente_lyapunov);
        fprintf(f4, "\n");
    //}
            
    fclose(archivo_posiciones);
    fclose(archivo_phi_psi);
    fclose(archivo_phi_dphi);
    fclose(f4);

     //Calculo el tiempo de ejecución
    //clock_t end=clock();
    //double tiempo= (double) (end-begin)/ CLOCKS_PER_SEC;
    //printf("Tiempo de compilacion=%lf \n", tiempo);

    return 0;
}