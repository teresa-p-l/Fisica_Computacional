import numpy as np
import matplotlib.pyplot as plt


pc = [ 0.069, 0.106, 0.144, 0.287]
pco1 = [0.026, 0.074, 0.106, 0.200]
pco2 = [0.026, 0.075, 0.104, 0.205]
pco3 = [0.058, 0.071, 0.097, 0.179]
pcofast = [0.056, 0.071, 0.094, 0.049]
pcomp = [0.034, 0.091, 0.145, 0.287]
joel = [0.053, 0.073, 0.112, 0.233]
joelo1 = [0.043, 0.064, 0.092, 0.182]
joelo2 = [0.046, 0.062, 0.093, 0.182]
joelo3 = [0.046 ,0.061, 0.09,0.176]
joelofast = [0.041, 0.062, 0.091, 0.130]
joelomp = [0.055, 0.073, 0.114, 0.234]



# Crear vector
N = [50, 100, 200, 500]


# Graficar
plt.figure(figsize=(10, 6))
plt.plot(N, pc, label=f'pc', color="#B266FF")
plt.plot(N, joel, label=f'joel', color="#9605BE")
plt.plot(N, pco1, label=f'pc o1', color="#FF0000")
plt.plot(N, joelo1, label=f'joel o1', color="#7F0404")
plt.plot(N, pco2, label=f'pc o2', color="#F8921D")
plt.plot(N, joelo2, label=f'joel o2', color="#FF4800")
plt.plot(N, pco3, label=f'pc o3', color="#49FE49")
plt.plot(N, joelo3, label=f'joel o3', color="#024802")
plt.plot(N, pcofast, label=f'pc ofast', color="#5E9EFF")
plt.plot(N, joelofast, label=f'joel ofast', color="#033BAA")
plt.plot(N, pcomp, label=f'pc omp', color="#FFFF00")
plt.plot(N, joelomp, label=f'joel omp', color="#858500")





plt.xticks(N, labels=[50, 100, 200, 500])
plt.title(f'Tiempos de ejecución con distintas optimizaciones en función del tiempo de ejecución')
plt.xlabel('Tiempo de ejecución (s)')
plt.ylabel('Tiempo de compilación')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Mostrar el gráfico por pantalla
plt.show()
