import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# === CONFIGURACIÓN ===
archivo1 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/perturbaciones/distancia_1.txt"  # cambia esto si el archivo tiene otro nombre
archivo3 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/perturbaciones/distancia_3.txt"  # cambia esto si el archivo tiene otro nombre
archivo5 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/perturbaciones/distancia_5.txt"  # cambia esto si el archivo tiene otro nombre
archivo10 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/perturbaciones/distancia_10.txt"  # cambia esto si el archivo tiene otro nombre
archivo15 = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/perturbaciones/distancia_15.txt"  # cambia esto si el archivo tiene otro nombre


dt = 0.01               # paso de integración usado en la simulación

# === LECTURA DEL ARCHIVO ===
delta1 = np.loadtxt(archivo1)          # lee los deltas: un delta por línea
delta3 = np.loadtxt(archivo3)
delta5 = np.loadtxt(archivo5)
delta10 = np.loadtxt(archivo10)
delta15 = np.loadtxt(archivo15)

t = np.arange(len(delta3)) * dt       # crea el vector de tiempos



def exponencial(x, a, b):
    return  0.000002*np.exp(a*x)

params1, cov = curve_fit(exponencial, t, delta1)
a1, b1 = params1
params3, cov = curve_fit(exponencial, t, delta3)
a3, b3 = params3
params5, cov = curve_fit(exponencial, t, delta5)
a5, b5 = params5
params10, cov = curve_fit(exponencial, t, delta10)
a10, b10 = params10
params15, cov = curve_fit(exponencial, t, delta15)
a15, b15 = params15




x1_fit = np.array([min(t), max(t)])
x3_fit = np.array([min(t), max(t)])
x5_fit = np.array([min(t), max(t)])
x10_fit = np.array([min(t), max(t)])
x15_fit = np.array([min(t), max(t)])

y1_fit = exponencial(t, a1, b1)
y3_fit = exponencial(t, a3, b3)
y5_fit = exponencial(t, a5, b5)
y10_fit = exponencial(t, a10, b10)
y15_fit = exponencial(t, a15, b15)



# === GRÁFICO 1: delta(t) ===
plt.figure(figsize=(6, 4))
transparencia = 0.5

plt.plot(t, delta1, label='E = 1.0', color='blue',alpha=transparencia)
plt.plot(t, y1_fit, color='blue', linestyle='-', linewidth=2.5, label= r'$\lambda$ = %f' % a1)
plt.plot(t, delta3, label='E = 3.0', color='orange',alpha=transparencia)
plt.plot(t, y3_fit, color='orange', linestyle='-', linewidth=2.5, label= r'$\lambda$ = %f' % a3)
plt.plot(t, delta5, label='E = 5.0', color='purple',alpha=transparencia)
plt.plot(t, y5_fit, color='purple', linestyle='-', linewidth=2.5, label= r'$\lambda$ = %f' % a5)
plt.plot(t, delta10, label='E = 10.0', color='red',alpha=transparencia)
plt.plot(t, y10_fit, color='red', linestyle='-', linewidth=2.5, label= r'$\lambda$ = %f' % a10)
#plt.plot(t, delta15, label='E = 15.0', color='yellow')


# AÑADIMOS LOS AJUSTES EXPONENCIALES




#plt.plot(x15_fit, y15_fit, color='yellow', linestyle='-', linewidth=2.5, label= r'$\lambda$ = %f' % a15)

plt.xlabel('Tiempo (t)')
plt.ylabel(r'$\delta(t)$')
plt.title('Separación entre trayectorias con $\Delta \phi$ = $\Delta \psi$ = 0.000001 ' )
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/plots/delta_t.png")
plt.show()

