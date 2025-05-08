import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

nombre_archivo = 'C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/momentos.txt'

# Leer el archivo y cargar las columnas
datos = np.loadtxt(nombre_archivo, skiprows=1)  # Cambia el delimitador si es necesario
x = datos[:, 1]  # Primera columna
y = datos[:, 2]  # Segunda columna

# Realizar el ajuste lineal
pendiente, intercepto, r_value, p_value, std_err = linregress(x, y)

    # Mostrar resultados
print(f"Pendiente: {pendiente}")
print(f"Intercepto: {intercepto}")
print(f"Coeficiente de correlación de Pearson (r): {r_value}")

    # Graficar los datos y la línea ajustada
plt.scatter(x, y, label='Datos', color='blue')
plt.plot(x, pendiente * x + intercepto, color='red', label='Ajuste lineal')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Ajuste Lineal')
plt.legend()
plt.grid()
plt.show()
