import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle, Rectangle
import os
E=10.0

# Función para leer el archivo de coordenadas
def leer_coordenadas(nombre_archivo):
    frames = []
    
    try:
        with open(nombre_archivo, 'r') as archivo:
            contenido = archivo.read()
            bloques = contenido.split('\n\n')
            
            for bloque in bloques:
                if bloque.strip():  # Ignorar bloques vacíos
                    lineas = bloque.strip().split('\n')
                    frame = []
                    
                    for linea in lineas:
                        if linea.strip():  # Ignorar líneas vacías
                            try:
                                x, y = map(float, linea.split(','))
                                frame.append((x, y))
                            except ValueError:
                                print(f"Error al parsear la línea: {linea}")
                    
                    if len(frame) == 2:  # Asegurarse de que hay exactamente 2 puntos
                        frames.append(frame)
        
        return frames
    except FileNotFoundError:
        print(f"El archivo {nombre_archivo} no se encontró.")
        return []
    except Exception as e:
        print(f"Error al leer el archivo: {e}")
        return []

# Función para crear la animación
def animar_pendulo_doble(frames):
    # Configuración de la figura
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title('Simulación de Péndulo Doble')
    
    # Inicializar elementos gráficos
    origen = (0, 0)
    radio_circulo = 0.1
    
    # Crear círculos
    circulo1 = Circle((0, 0), radio_circulo, fc='blue', zorder=10)
    circulo2 = Circle((0, 0), radio_circulo, fc='red', zorder=10)
    
    ax.add_patch(circulo1)
    ax.add_patch(circulo2)
    
    # Crear líneas
    linea1, = ax.plot([], [], 'k-', lw=2)
    linea2, = ax.plot([], [], 'k-', lw=2)
    
    # Función de inicialización para la animación
    def init():
        circulo1.center = (0, 0)
        circulo2.center = (0, 0)
        linea1.set_data([], [])
        linea2.set_data([], [])
        return circulo1, circulo2, linea1, linea2
    
    # Función de actualización para la animación
    def update(frame_num):
        if frame_num < len(frames):
            frame = frames[frame_num]
            
            # Actualizar posición de los círculos
            pos1 = frame[0]
            pos2 = frame[1]
            
            circulo1.center = pos1
            circulo2.center = pos2
            
            # Actualizar posición de las líneas
            linea1.set_data([origen[0], pos1[0]], [origen[1], pos1[1]])
            linea2.set_data([pos1[0], pos2[0]], [pos1[1], pos2[1]])
            
        return circulo1, circulo2, linea1, linea2
    
    # Crear la animación
    ani = animation.FuncAnimation(fig, update, frames=len(frames),
                                  init_func=init, blit=True, interval=10)

    return fig, ani

# Función principal

# Solicitar el nombre del archivo
nombre_archivo = f"C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/posicionesTXT/posiciones_E{E}.txt"
    
# Leer las coordenadas
frames = leer_coordenadas(nombre_archivo)
    

    
print(f"Se cargaron {len(frames)} frames de animación.")
    
# Crear la animación
fig, ani = animar_pendulo_doble(frames)
    
    # Guardar la animación como archivo GIF o mostrarla


    #ani.save(f"C:/Users/Teresa/Desktop/COMPU/Fisica_Computacional/Voluntario_PenduloDoble/animacion.mp4", writer='pillow', fps=20)
    #print(f"Animación guardada como {nombre_salida}.gif")

plt.show()
