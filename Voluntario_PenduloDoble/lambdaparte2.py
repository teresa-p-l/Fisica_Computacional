import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURACIÓN ===
archivo1 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_1.txt"  # cambia esto si el archivo tiene otro nombre
archivo3 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_3.txt"  # cambia esto si el archivo tiene otro nombre
archivo5 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_5.txt"  # cambia esto si el archivo tiene otro nombre
archivo10 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_10.txt"  # cambia esto si el archivo tiene otro nombre
archivo15 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_15.txt"  # cambia esto si el archivo tiene otro nombre
#archivo50 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/archivos/lambda_50.txt"  # cambia esto si el archivo tiene otro nombre


dt = 0.01               # paso de integración usado en la simulación

# === LECTURA DEL ARCHIVO ===
#delta1 = np.loadtxt(archivo1)          # lee los deltas: un delta por línea
delta3 = np.loadtxt(archivo3)
delta5 = np.loadtxt(archivo5)
delta10 = np.loadtxt(archivo10)
delta15 = np.loadtxt(archivo15)
#delta50 = np.loadtxt(archivo50)


t = np.arange(len(delta3)) * 10       # crea el vector de tiempos





# === GRÁFICO 1: delta(t) ===
plt.figure(figsize=(6, 4))
transparencia = 0.4

#plt.plot(t, delta1, label='E = 1.0', color='blue',alpha=transparencia)
plt.plot(t, delta3, label='E = 3.0', color='green',alpha=transparencia)
plt.plot(t, delta5, label='E = 5.0', color='purple',alpha=transparencia)
plt.plot(t, delta10, label='E = 10.0', color='red',alpha=transparencia)
plt.plot(t, delta15, label='E = 15.0', color='yellow',alpha=transparencia)
#plt.plot(t, delta50, label='E = 50.0', color='orange',alpha=transparencia)



plt.xlabel(r'Tiempo (t)')
plt.ylabel(r'$\lambda$')
plt.title('Separación entre trayectorias')
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/plots/delta_t_2cambiandoe.png")
plt.show()

