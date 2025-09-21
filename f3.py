import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import glob
import re
from scipy import ndimage
import warnings
from scipy import stats

# Ignorar advertencias para una ejecución más limpia
warnings.filterwarnings('ignore', category=np.RankWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Configurar directorios
input_dir = 'simulaciones'
output_dir = 'graficos'

# Crear directorios si no existen
os.makedirs(input_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

def parse_metadata(filename):
    """Extraer metadatos del archivo de simulación"""
    metadata = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    if ':' in line:
                        key, value = line[1:].split(':', 1)
                        metadata[key.strip()] = value.strip()
                else:
                    break
    except Exception as e:
        print(f"Error al leer el archivo {filename}: {str(e)}")
    return metadata

def load_simulation_data(filename):
    """Cargar datos de altura de la simulación"""
    try:
        print(f"Procesando archivo: {filename}")
        
        # Leer metadatos
        metadata = parse_metadata(filename)
        
        # Determinar tamaño de la rejilla desde metadatos
        grid_size = 50  # Valor por defecto
        for key in metadata:
            if 'Tamaño de la rejilla' in key:
                grid_match = re.search(r'(\d+)', metadata[key])
                if grid_match:
                    grid_size = int(grid_match.group(1))
                break
        
        # Leer datos de altura
        data = np.loadtxt(filename, comments='#')
        heights = np.zeros((grid_size, grid_size))
        
        for d in data:
            i, j, h = int(d[0])-1, int(d[1])-1, int(d[2])
            if i < grid_size and j < grid_size:
                heights[i, j] = h
        
        # Extraer información del nombre del archivo
        nombre_archivo = os.path.basename(filename)
        
        # Determinar el material
        if 'fe_' in nombre_archivo.lower() or 'fe' in nombre_archivo.lower():
            material = 'Fe'
        elif 'cr_' in nombre_archivo.lower() or 'cr' in nombre_archivo.lower():
            material = 'Cr'
        else:
            material = 'Desconocido'
        
        # Preparar resultados
        resultados = {
            'height': heights,
            'n_x': grid_size,
            'n_y': grid_size,
            'material': material,
            'filename': nombre_archivo
        }
        
        print(f"Procesado: {nombre_archivo} - Material: {material}")
        return resultados
        
    except Exception as e:
        print(f"Error al procesar {filename}: {str(e)}")
        return None

def generar_mapas_promedio_por_material(todos_datos):
    """Generar mapas promedio en 2D y 3D para cada material por separado"""
    if not todos_datos:
        print("No hay datos para generar mapas promedio")
        return
    
    # Separar datos por material
    datos_fe = [d for d in todos_datos if d['material'] == 'Fe']
    datos_cr = [d for d in todos_datos if d['material'] == 'Cr']
    
    # Procesar Fe
    if datos_fe:
        print("Generando mapas promedio para Fe...")
        generar_mapa_material(datos_fe, 'Fe')
    
    # Procesar Cr
    if datos_cr:
        print("Generando mapas promedio para Cr...")
        generar_mapa_material(datos_cr, 'Cr')

def generar_mapa_material(datos_material, nombre_material):
    """Generar mapas 2D y 3D para un material específico"""
    # Determinar el tamaño de la rejilla
    grid_size = datos_material[0]['n_x']
    
    # Crear matrices acumuladoras para el promedio
    height_sum = np.zeros((grid_size, grid_size))
    count = np.zeros((grid_size, grid_size))
    
    for datos in datos_material:
        height = datos['height']
        for i in range(grid_size):
            for j in range(grid_size):
                if i < height.shape[0] and j < height.shape[1]:
                    height_sum[i, j] += height[i, j]
                    count[i, j] += 1
    
    # Calcular el promedio (evitando división por cero)
    height_avg = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            if count[i, j] > 0:
                height_avg[i, j] = height_sum[i, j] / count[i, j]
    
    # Aplicar suavizado ligero
    height_smooth = ndimage.gaussian_filter(height_avg, sigma=0.5)
    
    # Generar mapa 2D con puntos (estilo gnuplot)
    plt.figure(figsize=(10, 8))
    
    # Crear coordenadas para el scatter plot
    x_coords = []
    y_coords = []
    z_values = []
    
    for i in range(grid_size):
        for j in range(grid_size):
            if height_smooth[i, j] > 0:  # Solo puntos con altura > 0
                x_coords.append(i)
                y_coords.append(j)
                z_values.append(height_smooth[i, j])
    
    # Crear scatter plot con puntos de tamaño variable según la altura
    scatter = plt.scatter(x_coords, y_coords, c=z_values, 
                         cmap='viridis', s=10, alpha=0.7,
                         marker='o', edgecolors='none')
    
    plt.colorbar(scatter, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa 2D Promedio de Altura - {nombre_material}\n{len(datos_material)} simulaciones', fontsize=14)
    plt.xlabel('Coordenada x (sitios)', fontsize=12)
    plt.ylabel('Coordenada y (sitios)', fontsize=12)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'mapa_2d_promedio_{nombre_material}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Generar mapa 3D
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    x = np.arange(grid_size)
    y = np.arange(grid_size)
    X, Y = np.meshgrid(x, y)
    
    # Usar un colormap apropiado
    cmap = plt.get_cmap('viridis')
    
    # Crear superficie 3D
    surf = ax.plot_surface(X, Y, height_smooth, 
                         cmap=cmap,
                         linewidth=0, 
                         antialiased=True, 
                         alpha=0.8,
                         rstride=1,
                         cstride=1)
    
    ax.set_zlim(0, np.max(height_smooth) * 1.1)
    ax.set_xlabel('Coordenada x (sitios)', fontsize=10, labelpad=10)
    ax.set_ylabel('Coordenada y (sitios)', fontsize=10, labelpad=10)
    ax.set_zlabel('Altura promedio (capas)', fontsize=10, labelpad=10)
    ax.set_title(f'Mapa 3D Promedio de Altura - {nombre_material}\n{len(datos_material)} simulaciones', fontsize=12, pad=20)
    ax.view_init(elev=30, azim=45)
    
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, pad=0.1)
    
    plt.savefig(os.path.join(output_dir, f'mapa_3d_promedio_{nombre_material}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Mapas para {nombre_material} generados correctamente")

def main():
    """Función principal"""
    print("Iniciando generación de mapas promedio...")
    
    # Verificar si hay archivos en el directorio de simulaciones
    archivos = glob.glob(os.path.join(input_dir, '*.dat'))
    if not archivos:
        print(f"No se encontraron archivos .dat en la carpeta '{input_dir}'")
        print("Por favor, coloca los archivos de simulación en la carpeta y ejecuta nuevamente")
        return
    
    print(f"Se encontraron {len(archivos)} archivos para procesar")
    
    # Procesar cada archivo
    todos_datos = []
    for archivo in archivos:
        resultados = load_simulation_data(archivo)
        if resultados is not None:
            todos_datos.append(resultados)
    
    if not todos_datos:
        print("No se pudieron procesar los archivos. Verifica el formato de los archivos .dat")
        return
    
    # Generar mapas promedio por material
    generar_mapas_promedio_por_material(todos_datos)
    
    print(f"Proceso completado. Los mapas promedio se han guardado en la carpeta '{output_dir}'")

if __name__ == "__main__":
    main()
