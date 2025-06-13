import numpy as np
import matplotlib.pyplot as plt




pc200 = [0.107, 0.105, 0.105, 0.105, 0.105]



joel200 = [0.115, 0.114, 0.114, 0.114, 0.114 ]

# Crear vector
N = [1, 4,8,12,16]


# Graficar
plt.figure(figsize=(10, 6))
plt.plot(N, pc200, label='pc tmax = 200s', color='red')
plt.plot(N, joel200, label='joel tmax = 200s', color='green')






plt.xticks(N, labels=N)
plt.title(f'Tiempo de compilación vs Número de Threads (para tmax = 200s)')
plt.xlabel('Threads')
plt.ylabel('Tiempo')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Mostrar el gráfico por pantalla
plt.show()