import numpy as np
import matplotlib.pyplot as plt

def maxwell_boltzmann_distribution(v, T, m, kB):
    """Calcula la distribución de velocidades de Maxwell-Boltzmann."""
    factor = (m / (kB * T)) 
    return factor * v * np.exp(-m * v**2 / (2 * kB * T))

# Parámetros
T = 290  # Temperatura en Kelvin
m = 1  # Masa de la partícula en kg (ejemplo: masa de un átomo de hidrógeno)
kB = 1  # Constante de Boltzmann en J/K

# Rango de velocidades
v = np.linspace(0, 100, 500)  # Velocidades en m/s

# Calcular la distribución
P_v = maxwell_boltzmann_distribution(v, T, m, kB)

# Leer datos del archivo
ruta_del_fichero = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/velocidades.txt"
datos = np.loadtxt(ruta_del_fichero)

# Graficar
plt.figure(figsize=(10, 6))

# Histograma de los datos
plt.hist(datos, bins=100, density=True, alpha=1, color='orange', label='Histograma de datos')

# Distribución de Maxwell-Boltzmann
plt.plot(v, P_v, label='Distribución de Maxwell-Boltzmann', color='blue')

# Configuración de la gráfica
plt.title('Distribución de Velocidades de Maxwell-Boltzmann y Histograma')
plt.xlabel('Velocidad (m/s)')
plt.ylabel('Densidad de probabilidad')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Mostrar la figura
plt.show()