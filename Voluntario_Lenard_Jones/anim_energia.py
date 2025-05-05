import numpy as np
import matplotlib.pyplot as plt

# Cargar los datos de energía
# Si el archivo tiene separadores por coma, cambia delimiter a ',' si es necesario
energia = np.loadtxt('C:/Users/Teresa/Desktop/COMPU/FISICA_COMPUTACIONAL/Fisica_Computacional/Voluntario_Lenard_Jones/energia.txt')

# Separar columnas
ekin = energia[:, 0]  # Energía cinética
epot = energia[:, 1]  # Energía potencial
etot = energia[:, 2]  # Energía total

# Crear vector de tiempo / frames
frames = np.arange(len(energia))

# Graficar
plt.figure(figsize=(10, 6))
plt.plot(frames, ekin, label='Energía Cinética', color='r')
plt.plot(frames, epot, label='Energía Potencial', color='g')
plt.plot(frames, etot, label='Energía Total', color='b')

plt.title('Energías del Sistema')
plt.xlabel('Frame')
plt.ylabel('Energía')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Mostrar el gráfico por pantalla
plt.show()
