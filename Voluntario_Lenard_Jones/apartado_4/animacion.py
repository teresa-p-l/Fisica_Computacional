import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


data = np.loadtxt('C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_4/datos_simulacion.txt')
N, L, h, T = int(data[0]), data[1], data[2], data[3]

# Parámetros a adaptar
filename    = "C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_4/posiciones.txt"   # nombre de tu fichero              
interval_ms = 2                # tiempo entre frames en ms
output_mp4  = "C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_4/simulacion.mp4"   # nombre del fichero de salida

# --- 1) Leer y estructurar datos ---
# Cada bloque de N líneas, dos columnas x,y
data = np.loadtxt(filename)
n_frames = data.shape[0] // N

# reshaped: (n_frames, N, 2)
frames = data.reshape((n_frames, N, 2))

# --- 2) Preparar figura ---
fig, ax = plt.subplots()
scat = ax.scatter([], [], s=50)

# Ajusta límites en función de tus datos
ax.set_xlim( 0, L)
ax.set_ylim(0, L)
ax.set_xlabel('x'); ax.set_ylabel('y')

# Función de inicialización
def init():
    scat.set_offsets(np.empty((0,2)))
    return scat,

# Función que actualiza cada frame
def update(frame_idx):
    xy = frames[frame_idx]
    scat.set_offsets(xy)
    ax.set_title(f'Tiempo: frame {frame_idx+1}/{n_frames}')
    return scat,

# --- 3) Construir animación ---
anim = animation.FuncAnimation(fig, update,
                               frames=range(n_frames),
                               init_func=init,
                               interval=interval_ms,
                               blit=True)

# Mostrar por pantalla
plt.show()

# Guardar como vídeo mp4 (requiere ffmpeg instalado)
#anim.save(output_mp4, writer='ffmpeg', dpi=150)
#print(f'Vídeo guardado en {output_mp4}')










# Cargar los datos de energía
# Si el archivo tiene separadores por coma, cambia delimiter a ',' si es necesario
energia = np.loadtxt('C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_4/energia.txt')

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













# Hacemos la ecuación de maxwell

import numpy as np
import matplotlib.pyplot as plt

def maxwell_boltzmann_distribution(v, T, m, kB):
    """Calcula la distribución de velocidades de Maxwell-Boltzmann."""
    factor = (m / (kB * T)) 
    return factor * v * np.exp(-m * v**2 / (2 * kB * T))

# Parámetros
m = 1  # Masa de la partícula en kg (ejemplo: masa de un átomo de hidrógeno)
kB = 1  # Constante de Boltzmann en J/K



# Leer datos del archivo
ruta_del_fichero = "c:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_4/velocidades.txt"
datos = np.loadtxt(ruta_del_fichero)

# Rango de velocidades
v = np.linspace(0, np.max(datos)*1.5, 200)  # Velocidades en m/s

# Calcular la distribución
P_v = maxwell_boltzmann_distribution(v, T, m, kB)

# Graficar
plt.figure(figsize=(10, 6))

# Histograma de los datos
plt.hist(datos, bins=30, density=True, alpha=1, color='orange', label='Histograma de datos')

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





