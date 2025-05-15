import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


data = np.loadtxt('C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/datos_simulacion.txt')
N, L, h, T = int(data[0]), data[1], data[2], data[3]

# Parámetros a adaptar
filename    = "C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/posiciones.txt"   # nombre de tu fichero              
interval_ms = 2                # tiempo entre frames en ms
output_mp4  = "C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/simulacion.mp4"   # nombre del fichero de salida

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
energia = np.loadtxt('C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/energia.txt')

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





















# APARTADO 7: VIDEO DE LAS FLUCTUACIONES


#CAMBIA LA PARTICULA QUE FIJAS

particula=0

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

filename = 'C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/desplazamiento.txt'

# Función para leer las fluctuaciones desde el archivo
def leer_fluctuaciones(filename):

    fluctuaciones = []
    with open(filename, 'r') as file:
        particulas_actuales = []  # Almacena las fluctuaciones de las partículas en un paso temporal
        for linea in file:
            linea = linea.strip()
            if linea:  # Si la línea no está vacía
                datos = list(map(float, linea.split()))  # Convertir los valores a flotantes
                particulas_actuales.append(datos)
            else:  # Si encontramos un salto de línea, procesamos el paso temporal
                if particulas_actuales:
                    # Transponer para obtener la evolución temporal de cada partícula
                    if not fluctuaciones:
                        fluctuaciones = [[] for _ in range(len(particulas_actuales))]
                    for i, particula in enumerate(particulas_actuales):
                        fluctuaciones[i].extend(particula)
                    particulas_actuales = []  # Reiniciar para el siguiente paso temporal

        # Procesar el último paso si no está vacío
        if particulas_actuales:
            if not fluctuaciones:
                fluctuaciones = [[] for _ in range(len(particulas_actuales))]
            for i, particula in enumerate(particulas_actuales):
                fluctuaciones[i].extend(particula)

    return fluctuaciones

# Función para crear el video directamente
def crear_video(fluctuaciones, video_path, intervalo=500):

    writer = FFMpegWriter(fps=1000 // intervalo)
    fig, ax = plt.subplots(figsize=(8, 6))

    with writer.saving(fig, video_path, dpi=100):
        for particula_idx, fluctuacion in enumerate(fluctuaciones):
            pasos = range(1, len(fluctuacion) + 1)
            ax.clear()
            ax.plot(pasos, fluctuacion, color='blue', linewidth=0.8, label=f'Partícula {particula_idx + 1}')
            ax.set_xlabel('Paso temporal')
            ax.set_ylabel('Fluctuación')
            ax.set_title(f'Fluctuación de la Partícula {particula_idx + 1} con respecto a la partícula {particula}')
            ax.legend()
            ax.grid(True)
            writer.grab_frame()

    print(f"Video guardado como '{video_path}'.")

# Ejemplo de uso
if __name__ == "__main__":
    archivo = 'C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/desplazamiento.txt'  # Nombre del archivo de fluctuaciones
    video_path = f"C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_Lenard_Jones/apartado_7/fluctuaciones_{particula}.mp4"

    # Leer las fluctuaciones
    fluctuaciones = leer_fluctuaciones(archivo)

    # Crear el video directamente
    crear_video(fluctuaciones, video_path, intervalo=500)