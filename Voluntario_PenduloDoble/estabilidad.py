import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import linregress,pearsonr
from scipy.optimize import curve_fit



def convert_float(value):
    # Convert bytes to string if necessary
    if isinstance(value, bytes):
        value = value.decode('utf-8')
    return float(value.replace(',', '.'))

def load_data_from_file(filename):
    # Define converters for each column to handle commas as decimal separators
    converters = {0: convert_float, 1: convert_float}#, 2: convert_float, 3: convert_float}
    
    # Load the data from the file using the converters and specifying the tab delimiter
    data = np.loadtxt(filename, delimiter=',', converters=converters)
    
    # Extract each column from the data
    x = data[:, 0]
    y = data[:, 1]
    #errx = data[:, 2]
    #erry = data[:, 3]
    
    return x, y#, errx, erry

# Cambiar nombre del archivo de datos para cambiar los datos recogidos
x, y = load_data_from_file("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/poincareTXT/poincare_phi_psi_Estabilidad1.0.txt")
x1, y1 = load_data_from_file("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/poincareTXT/poincare_phi_psi_Estabilidad2.0.txt")
x2, y2 = load_data_from_file("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/poincareTXT/poincare_phi_psi_Estabilidad3.0.txt")
x3, y3  = load_data_from_file("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/poincareTXT/poincare_phi_psi_Estabilidad4.0.txt")
x4, y4 = load_data_from_file("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/poincareTXT/poincare_phi_psi_Estabilidad5.0.txt")
#x5, y5 = load_data_from_file("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/poincareTXT/poincare_phi_psi_E20.0.txt")

#def line_fit(x, slope, intercept):
  #  return slope * x + intercept

#Valores de x para la línea ajustada
line_x = np.array([min(x), max(x)])
line_x1 = np.array([min(x1), max(x1)])
line_x2 = np.array([min(x2), max(x2)])
line_x3 = np.array([min(x3), max(x3)])
line_x4 = np.array([min(x4), max(x4)])
#line_x5 = np.array([min(x5), max(x5)])
# Valores de y usando la función de ajuste
line_y = np.array([min(y), max(y)])
line_y1 = np.array([min(y1), max(y1)])
line_y2 = np.array([min(y2), max(y2)])
line_y3 = np.array([min(y3), max(y3)])
line_y4 = np.array([min(y4), max(y4)])
#line_y5 = np.array([min(y5), max(y5)])



#Cálculo del coef pearson
pearson_coef, _ = pearsonr(x, y)
print (pearson_coef)

num_bins = 50

plt.figure(figsize=(8, 6))



#plt.plot(x3, y3, color='orange', label='$\phi$ = 0.4 rad', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x2, y2, color='purple', label='$\phi$ = 0.25 rad', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x, y, color='blue', label='$\phi$ = 0.22 rad', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x1, y1, color='red', label='$\phi$ = 0.2 rad', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x4, y4, color='green', label='$\phi$ = 0.1 rad', alpha=0.5)  # Dibuja los puntos de datos con transparencia
#plt.plot(x5, y5, color='yellow', label='E=20', alpha=0.5)  # Dibuja los puntos de datos con transparencia



plt.xlabel('$\phi$ (rad)')
plt.ylabel('$\psi$ (rad)')
plt.title('Mapa de Poincaré variando $\phi$ y $\psi$=0.2 rad y E=1.0')




plt.legend()
plt.grid(True)  #Poner cuadrilla
plt.savefig("C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/plots/estabilidad2")  #Guardar plot en ruta
plt.show()

#se_intercept = std_err * np.sqrt(np.sum(x**2) / len(x))
#print(f"La ecuación de la línea es y = ({slope:.9f} ± {std_err:.9f})x + ({intercept:.9f} ± {se_intercept:.9f})")