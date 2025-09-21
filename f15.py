import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re
from scipy import ndimage
import warnings
from matplotlib.colors import LinearSegmentedColormap

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
        
        return resultados
        
    except Exception as e:
        print(f"Error al procesar {filename}: {str(e)}")
        return None

def generar_mapas_promedio_2d(datos_material, nombre_material):
    """Generar múltiples estilos de mapas 2D para un material específico"""
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
    
    # Crear directorio para el material si no existe
    material_dir = os.path.join(output_dir, nombre_material)
    os.makedirs(material_dir, exist_ok=True)
    
    # 1. Mapa de puntos (estilo gnuplot)
    generar_mapa_puntos(height_smooth, grid_size, nombre_material, material_dir)
    
    # 2. Mapa de calor (heatmap)
    generar_mapa_calor(height_smooth, grid_size, nombre_material, material_dir)
    
    # 3. Mapa de contorno (contour)
    generar_mapa_contorno(height_smooth, grid_size, nombre_material, material_dir)
    
    # 4. Mapa de superficie 3D (proyección 2D)
    generar_mapa_superficie(height_smooth, grid_size, nombre_material, material_dir)
    
    # 5. Mapa granulado (con textura)
    generar_mapa_granulado(height_smooth, grid_size, nombre_material, material_dir)
    
    # 6. Mapa con interpolación (suavizado)
    generar_mapa_interpolado(height_smooth, grid_size, nombre_material, material_dir)

def generar_mapa_puntos(height_data, grid_size, material, output_dir):
    """Generar mapa de puntos estilo gnuplot"""
    plt.figure(figsize=(10, 8))
    
    # Crear coordenadas para el scatter plot
    x_coords = []
    y_coords = []
    z_values = []
    
    for i in range(grid_size):
        for j in range(grid_size):
            if height_data[i, j] > 0:  # Solo puntos con altura > 0
                x_coords.append(i)
                y_coords.append(j)
                z_values.append(height_data[i, j])
    
    # Crear scatter plot con puntos de tamaño variable según la altura
    scatter = plt.scatter(x_coords, y_coords, c=z_values, 
                         cmap='viridis', s=15, alpha=0.8,
                         marker='o', edgecolors='none')
    
    plt.colorbar(scatter, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa 2D de Puntos - {material}\n{len(x_coords)} puntos', fontsize=14)
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
    plt.savefig(os.path.join(output_dir, f'mapa_puntos_{material}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generar_mapa_calor(height_data, grid_size, material, output_dir):
    """Generar mapa de calor"""
    plt.figure(figsize=(10, 8))
    
    # Crear mapa de calor
    im = plt.imshow(height_data, 
                   cmap='hot',
                   interpolation='nearest',
                   origin='lower',
                   aspect='equal')
    
    plt.colorbar(im, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa de Calor - {material}', fontsize=14)
    plt.xlabel('Coordenada x (sitios)', fontsize=12)
    plt.ylabel('Coordenada y (sitios)', fontsize=12)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'mapa_calor_{material}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generar_mapa_contorno(height_data, grid_size, material, output_dir):
    """Generar mapa de contorno"""
    plt.figure(figsize=(10, 8))
    
    # Crear coordenadas
    x = np.arange(grid_size)
    y = np.arange(grid_size)
    X, Y = np.meshgrid(x, y)
    
    # Crear mapa de contorno
    contour = plt.contourf(X, Y, height_data, 20, cmap='plasma')
    
    plt.colorbar(contour, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa de Contorno - {material}', fontsize=14)
    plt.xlabel('Coordenada x (sitios)', fontsize=12)
    plt.ylabel('Coordenada y (sitios)', fontsize=12)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'mapa_contorno_{material}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generar_mapa_superficie(height_data, grid_size, material, output_dir):
    """Generar mapa de superficie (proyección 2D)"""
    plt.figure(figsize=(10, 8))
    
    # Crear coordenadas
    x = np.arange(grid_size)
    y = np.arange(grid_size)
    X, Y = np.meshgrid(x, y)
    
    # Crear mapa de superficie con sombreado
    im = plt.pcolormesh(X, Y, height_data, cmap='terrain', shading='auto')
    
    plt.colorbar(im, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa de Superficie - {material}', fontsize=14)
    plt.xlabel('Coordenada x (sitios)', fontsize=12)
    plt.ylabel('Coordenada y (sitios)', fontsize=12)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'mapa_superficie_{material}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generar_mapa_granulado(height_data, grid_size, material, output_dir):
    """Generar mapa con textura granulado"""
    plt.figure(figsize=(10, 8))
    
    # Añadir ruido para efecto granulado
    np.random.seed(42)  # Para reproducibilidad
    noise = np.random.normal(0, 0.1, height_data.shape)
    height_noisy = np.clip(height_data + noise, 0, None)
    
    # Crear mapa con textura
    im = plt.imshow(height_noisy, 
                   cmap='viridis',
                   interpolation='nearest',
                   origin='lower',
                   aspect='equal')
    
    plt.colorbar(im, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa Granulado - {material}', fontsize=14)
    plt.xlabel('Coordenada x (sitios)', fontsize=12)
    plt.ylabel('Coordenada y (sitios)', fontsize=12)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'mapa_granulado_{material}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generar_mapa_interpolado(height_data, grid_size, material, output_dir):
    """Generar mapa con interpolación suavizada"""
    plt.figure(figsize=(10, 8))
    
    # Crear mapa con interpolación
    im = plt.imshow(height_data, 
                   cmap='coolwarm',
                   interpolation='bicubic',
                   origin='lower',
                   aspect='equal')
    
    plt.colorbar(im, label='Altura promedio (capas)', shrink=0.8)
    
    plt.title(f'Mapa Interpolado - {material}', fontsize=14)
    plt.xlabel('Coordenada x (sitios)', fontsize=12)
    plt.ylabel('Coordenada y (sitios)', fontsize=12)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'mapa_interpolado_{material}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Función principal"""
    print("Iniciando generación de mapas promedio 2D...")
    
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
            print(f"Procesado: {resultados['filename']} - Material: {resultados['material']}")
    
    if not todos_datos:
        print("No se pudieron procesar los archivos. Verifica el formato de los archivos .dat")
        return
    
    # Separar datos por material
    datos_fe = [d for d in todos_datos if d['material'] == 'Fe']
    datos_cr = [d for d in todos_datos if d['material'] == 'Cr']
    
    # Generar mapas para Fe
    if datos_fe:
        print(f"Generando mapas para Fe ({len(datos_fe)} archivos)...")
        generar_mapas_promedio_2d(datos_fe, 'Fe')
    
    # Generar mapas para Cr
    if datos_cr:
        print(f"Generando mapas para Cr ({len(datos_cr)} archivos)...")
        generar_mapas_promedio_2d(datos_cr, 'Cr')
    
    print(f"Proceso completado. Los mapas se han guardado en la carpeta '{output_dir}'")

if __name__ == "__main__":
    main()
