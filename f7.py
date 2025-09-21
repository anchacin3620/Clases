import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import glob
import re
from scipy import ndimage
import warnings
from scipy import stats
import csv

# Ignorar advertencias para una ejecución más limpia
warnings.filterwarnings('ignore', category=np.RankWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Configurar directorios
input_dir = 'simulaciones'
output_dir = 'graficos'

# Crear directorios si no existen
os.makedirs(input_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Definir umbrales para identificar etapas
UMBRALES_ETAPAS = {
    'nucleacion': 1.5,
    'coalescencia': (1.5, 5.0),
    'crecimiento': 5.0
}

def identificar_etapa(datos):
    """Identifica la etapa de crecimiento considerando espesor, cobertura y densidad de islas."""
    espesor = datos['espesor_promedio']
    cobertura = datos.get('cobertura_fraccional', 0)
    densidad_islas = datos.get('densidad_islas', 0)
    tamano_promedio = datos.get('tamano_promedio_islas', 0)
    
    # Etapa de nucleación: baja cobertura, alta densidad de islas pequeñas
    if espesor < UMBRALES_ETAPAS['nucleacion'] and cobertura < 20 and densidad_islas > 0.005 and tamano_promedio < 10:
        return 'nucleacion'
    # Etapa de coalescencia: cobertura media, densidad de islas disminuyendo
    elif UMBRALES_ETAPAS['coalescencia'][0] <= espesor <= UMBRALES_ETAPAS['coalescencia'][1] and cobertura < 70:
        return 'coalescencia'
    # Etapa de crecimiento: alta cobertura, baja densidad de islas
    else:
        return 'crecimiento'

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
        
        # Calcular estadísticas básicas
        rugosidad = np.std(heights)
        espesor_promedio = np.mean(heights)
        max_altura = np.max(heights)
        min_altura = np.min(heights)
        rugosidad_rms = np.sqrt(np.mean(heights**2))
        
        # Calcular asimetría y curtosis
        altura_aplanada = heights.flatten()
        asimetria = stats.skew(altura_aplanada)
        curtosis = stats.kurtosis(altura_aplanada)
        
        # Extraer información del nombre del archivo
        nombre_archivo = os.path.basename(filename)
        
        # Determinar el material
        if 'fe_' in nombre_archivo.lower() or 'fe' in nombre_archivo.lower():
            material = 'Fe'
        elif 'cr_' in nombre_archivo.lower() or 'cr' in nombre_archivo.lower():
            material = 'Cr'
        else:
            material = 'Desconocido'
        
        # Extraer parámetros de simulación
        energia_difusion = 0.0
        tasa_eventos = 0.0
        energia_promedio = 0.0
        
        # Buscar parámetros en el nombre del archivo
        match_ed = re.search(r'ed([\d.]+)', nombre_archivo)
        if match_ed:
            energia_difusion = float(match_ed.group(1))
        
        match_te = re.search(r'te([\d.]+)', nombre_archivo)
        if match_te:
            tasa_eventos = float(match_te.group(1))
        
        match_ep = re.search(r'ep([\d.]+)', nombre_archivo)
        if match_ep:
            energia_promedio = float(match_ep.group(1))
        
        # Calcular métricas adicionales
        estructura = ndimage.generate_binary_structure(2, 2)
        mascara = heights > 2  # Solo material depositado, excluyendo sustrato
        islas, num_islas = ndimage.label(mascara, structure=estructura)
        densidad_islas = num_islas / (grid_size * grid_size) if grid_size * grid_size > 0 else 0
        cobertura_fraccional = np.sum(heights > 2) / (grid_size * grid_size) * 100 if grid_size * grid_size > 0 else 0
        
        # Calcular cobertura por umbrales
        umbrales = [10, 20, 30, 60]
        coberturas_umbrales = {}
        for umbral in umbrales:
            areas_depositadas = heights > umbral
            cobertura = np.sum(areas_depositadas) / (grid_size * grid_size) * 100 if grid_size * grid_size > 0 else 0
            coberturas_umbrales[f'cobertura_{umbral}'] = cobertura
        
        # Identificar etapa de crecimiento
        etapa = identificar_etapa({
            'espesor_promedio': espesor_promedio,
            'cobertura_fraccional': cobertura_fraccional,
            'densidad_islas': densidad_islas,
            'tamano_promedio_islas': 0  # No calculado en esta versión simplificada
        })
        
        # Preparar resultados
        resultados = {
            'height': heights,
            'rugosidad': rugosidad, 
            'rugosidad_rms': rugosidad_rms,
            'asimetria': asimetria,
            'curtosis': curtosis,
            'espesor_promedio': espesor_promedio,
            'altura_maxima': max_altura, 
            'altura_minima': min_altura,
            'n_x': grid_size,
            'n_y': grid_size,
            'material': material,
            'energia_difusion': energia_difusion,
            'tasa_eventos': tasa_eventos,
            'energia_promedio': energia_promedio,
            'num_islas': num_islas,
            'densidad_islas': densidad_islas,
            'cobertura_fraccional': cobertura_fraccional,
            'coberturas_umbrales': coberturas_umbrales,
            'etapa': etapa,
            'filename': nombre_archivo
        }
        
        print(f"Procesado: {nombre_archivo} - Material: {material} - Rugosidad: {rugosidad:.4f}")
        return resultados
        
    except Exception as e:
        print(f"Error al procesar {filename}: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

def create_visualizations(filename, resultados):
    """Crear visualizaciones 2D, 3D e histogramas para un archivo individual"""
    try:
        height = resultados['height']
        material = resultados.get('material', 'Desconocido')
        energia_dif = resultados.get('energia_difusion', 'N/A')
        rugosidad = resultados.get('rugosidad', 0)
        rugosidad_rms = resultados.get('rugosidad_rms', 0)
        asimetria = resultados.get('asimetria', 0)
        curtosis = resultados.get('curtosis', 0)
        
        nombre_base = os.path.basename(filename).replace(".dat", "")
        
        # Visualización 2D - MEJORADA
        plt.figure(figsize=(12, 10))
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['axes.titlesize'] = 16
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        
        # Usar un colormap de alta calidad
        cmap = plt.get_cmap('viridis')
        
        # Aplicar suavizado a los datos para mejor visualización
        height_smooth = ndimage.gaussian_filter(height, sigma=1.0)
        
        # Crear el gráfico con interpolación para mejor calidad
        im = plt.imshow(height_smooth, 
                       cmap=cmap,
                       interpolation='bilinear',  # Interpolación bilineal para suavizar
                       vmin=0, 
                       vmax=resultados['altura_maxima'],
                       origin='lower',
                       aspect='equal')  # Mantener relación de aspecto 1:1
        
        # Mejorar la barra de color
        cbar = plt.colorbar(im, label='Altura (capas)', shrink=0.8, pad=0.03)
        cbar.ax.tick_params(labelsize=12)
        
        plt.title(f'Mapa 2D de Altura - {material}\nEnergía de difusión: {energia_dif} eV', fontsize=16, pad=20)
        plt.xlabel('Coordenada x (sitios)', fontsize=14, labelpad=10)
        plt.ylabel('Coordenada y (sitios)', fontsize=14, labelpad=10)
        
        # Ajustar los ticks para mejor visualización
        n_ticks = 6
        xticks = np.linspace(0, resultados['n_x']-1, n_ticks, dtype=int)
        yticks = np.linspace(0, resultados['n_y']-1, n_ticks, dtype=int)
        plt.xticks(xticks, xticks)
        plt.yticks(yticks, yticks)
        
        # Añadir grid para mejor referencia
        plt.grid(False)  # Desactivar grid para no interferir con la visualización
        
        # Añadir texto con información adicional
        textstr = f'Rugosidad: {rugosidad:.3f}\nEspesor prom: {resultados["espesor_promedio"]:.1f}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'mapa_2d_{nombre_base}.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Visualización 3D
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        x = np.arange(resultados['n_x'])
        y = np.arange(resultados['n_y'])
        X, Y = np.meshgrid(x, y)
        
        # Suavizar los datos para la visualización 3D
        height_smooth_3d = ndimage.gaussian_filter(height, sigma=1.0)
        
        surf = ax.plot_surface(X, Y, height_smooth_3d, 
                             cmap=cmap,
                             linewidth=0, 
                             antialiased=True, 
                             alpha=0.8,
                             rstride=1,  # Reducir el salto para mayor detalle
                             cstride=1)  # Reducir el salto para mayor detalle
        
        ax.set_zlim(0, resultados['altura_maxima'] * 1.1)
        ax.set_xlabel('Coordenada x (sitios)', fontsize=12, labelpad=10)
        ax.set_ylabel('Coordenada y (sitios)', fontsize=12, labelpad=10)
        ax.set_zlabel('Altura (capas)', fontsize=12, labelpad=10)
        ax.set_title(f'Mapa 3D de Altura - {material}\nEnergía de difusión: {energia_dif} eV', fontsize=14, pad=20)
        ax.view_init(elev=30, azim=45)
        
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, pad=0.1)
        
        plt.savefig(os.path.join(output_dir, f'mapa_3d_{nombre_base}.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Histograma de alturas
        plt.figure(figsize=(10, 8))
        plt.hist(height.flatten(), bins=30, alpha=0.7, edgecolor='black', density=True)
        plt.title(f'Distribución de Alturas - {material}\nEnergía de difusión: {energia_dif} eV', fontsize=16)
        plt.xlabel('Altura (capas atómicas)', fontsize=14)
        plt.ylabel('Densidad de probabilidad', fontsize=14)
        plt.grid(True, alpha=0.3)
        
        # Añadir estadísticas al gráfico
        textstr = f'Media: {resultados["espesor_promedio"]:.2f}\nDesviación: {rugosidad:.2f}\nAsimetría: {asimetria:.2f}\nCurtosis: {curtosis:.2f}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        plt.text(0.75, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        
        plt.savefig(os.path.join(output_dir, f'histograma_{nombre_base}.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Gráficos generados para: {nombre_base}")
        
    except Exception as e:
        print(f"Error al generar gráficos para {filename}: {str(e)}")
        plt.close()

def generar_graficos_comparativos(todos_datos):
    """Generar gráficos comparativos entre todos los datos"""
    if not todos_datos:
        print("No hay datos para generar gráficos comparativos")
        return
    
    # Separar datos por material
    datos_fe = [d for d in todos_datos if d['material'] == 'Fe']
    datos_cr = [d for d in todos_datos if d['material'] == 'Cr']
    
    # Crear gráfico comparativo de rugosidad
    plt.figure(figsize=(12, 8))
    plt.rcParams['font.size'] = 12
    
    if datos_fe:
        rugosidades_fe = [d['rugosidad'] for d in datos_fe]
        plt.hist(rugosidades_fe, alpha=0.5, label='Fe', bins=15, density=True)
    
    if datos_cr:
        rugosidades_cr = [d['rugosidad'] for d in datos_cr]
        plt.hist(rugosidades_cr, alpha=0.5, label='Cr', bins=15, density=True)
    
    plt.xlabel('Rugosidad (capas)', fontsize=14)
    plt.ylabel('Densidad de probabilidad', fontsize=14)
    plt.title('Comparación de Rugosidad: Fe vs Cr', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'comparacion_rugosidad.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Crear gráfico comparativo de espesor
    plt.figure(figsize=(12, 8))
    
    if datos_fe:
        espesores_fe = [d['espesor_promedio'] for d in datos_fe]
        plt.hist(espesores_fe, alpha=0.5, label='Fe', bins=15, density=True)
    
    if datos_cr:
        espesores_cr = [d['espesor_promedio'] for d in datos_cr]
        plt.hist(espesores_cr, alpha=0.5, label='Cr', bins=15, density=True)
    
    plt.xlabel('Espesor promedio (capas)', fontsize=14)
    plt.ylabel('Densidad de probabilidad', fontsize=14)
    plt.title('Comparación de Espesor: Fe vs Cr', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'comparacion_espesor.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Crear gráfico de dispersión: Rugosidad vs Energía de difusión
    plt.figure(figsize=(12, 8))
    
    if datos_fe:
        energias_fe = [d['energia_difusion'] for d in datos_fe]
        rugosidades_fe = [d['rugosidad'] for d in datos_fe]
        plt.scatter(energias_fe, rugosidades_fe, alpha=0.7, label='Fe', s=80)
    
    if datos_cr:
        energias_cr = [d['energia_difusion'] for d in datos_cr]
        rugosidades_cr = [d['rugosidad'] for d in datos_cr]
        plt.scatter(energias_cr, rugosidades_cr, alpha=0.7, label='Cr', s=80)
    
    plt.xlabel('Energía de difusión (eV)', fontsize=14)
    plt.ylabel('Rugosidad (capas)', fontsize=14)
    plt.title('Rugosidad vs Energía de difusión', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'rugosidad_vs_energia.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Gráficos comparativos generados")

def generar_graficos_promedio(todos_datos):
    """Generar gráficos promedio de todos los datos"""
    if not todos_datos:
        print("No hay datos para generar gráficos promedio")
        return
    
    # Determinar el tamaño de la rejilla
    grid_size = todos_datos[0]['n_x']
    
    # Crear matriz promedio
    height_sum = np.zeros((grid_size, grid_size))
    count = np.zeros((grid_size, grid_size))
    
    for datos in todos_datos:
        height = datos['height']
        for i in range(grid_size):
            for j in range(grid_size):
                if i < height.shape[0] and j < height.shape[1]:
                    height_sum[i, j] += height[i, j]
                    count[i, j] += 1
    
    # Calcular promedio
    height_avg = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            if count[i, j] > 0:
                height_avg[i, j] = height_sum[i, j] / count[i, j]
    
    # Aplicar suavizado
    height_smooth = ndimage.gaussian_filter(height_avg, sigma=1.0)
    
    # Gráfico 2D promedio
    plt.figure(figsize=(12, 10))
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    
    cmap = plt.get_cmap('viridis')
    
    im = plt.imshow(height_smooth, 
                   cmap=cmap,
                   interpolation='bilinear',
                   vmin=0, 
                   vmax=np.max(height_avg),
                   origin='lower',
                   aspect='equal')
    
    cbar = plt.colorbar(im, label='Altura promedio (capas)', shrink=0.8, pad=0.03)
    cbar.ax.tick_params(labelsize=12)
    
    plt.title(f'Mapa 2D Promedio de Altura\n{len(todos_datos)} simulaciones', fontsize=16, pad=20)
    plt.xlabel('Coordenada x (sitios)', fontsize=14, labelpad=10)
    plt.ylabel('Coordenada y (sitios)', fontsize=14, labelpad=10)
    
    # Ajustar los ticks
    n_ticks = 6
    xticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    yticks = np.linspace(0, grid_size-1, n_ticks, dtype=int)
    plt.xticks(xticks, xticks)
    plt.yticks(yticks, yticks)
    
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mapa_2d_promedio.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gráfico 3D promedio
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    x = np.arange(grid_size)
    y = np.arange(grid_size)
    X, Y = np.meshgrid(x, y)
    
    height_smooth_3d = ndimage.gaussian_filter(height_avg, sigma=1.0)
    
    surf = ax.plot_surface(X, Y, height_smooth_3d, 
                         cmap=cmap,
                         linewidth=0, 
                         antialiased=True, 
                         alpha=0.8,
                         rstride=1,
                         cstride=1)
    
    ax.set_zlim(0, np.max(height_avg) * 1.1)
    ax.set_xlabel('Coordenada x (sitios)', fontsize=12, labelpad=10)
    ax.set_ylabel('Coordenada y (sitios)', fontsize=12, labelpad=10)
    ax.set_zlabel('Altura promedio (capas)', fontsize=12, labelpad=10)
    ax.set_title(f'Mapa 3D Promedio de Altura\n{len(todos_datos)} simulaciones', fontsize=14, pad=20)
    ax.view_init(elev=30, azim=45)
    
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, pad=0.1)
    
    plt.savefig(os.path.join(output_dir, 'mapa_3d_promedio.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Gráficos promedio generados")

def guardar_resultados_csv(todos_datos):
    """Guardar todos los resultados en un archivo CSV usando el módulo csv"""
    if not todos_datos:
        print("No hay datos para guardar en CSV")
        return
    
    # Preparar datos para CSV
    campos = ['archivo', 'material', 'rugosidad', 'rugosidad_rms', 'asimetria', 'curtosis',
              'espesor_promedio', 'altura_maxima', 'altura_minima', 'num_islas',
              'densidad_islas', 'cobertura_fraccional', 'energia_difusion',
              'tasa_eventos', 'energia_promedio', 'etapa']
    
    # Añadir campos para coberturas por umbral
    for umbral in [10, 20, 30, 60]:
        campos.append(f'cobertura_{umbral}')
    
    csv_path = os.path.join(output_dir, 'resultados_analisis.csv')
    
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=campos)
        writer.writeheader()
        
        for datos in todos_datos:
            fila = {
                'archivo': datos['filename'],
                'material': datos['material'],
                'rugosidad': datos['rugosidad'],
                'rugosidad_rms': datos['rugosidad_rms'],
                'asimetria': datos['asimetria'],
                'curtosis': datos['curtosis'],
                'espesor_promedio': datos['espesor_promedio'],
                'altura_maxima': datos['altura_maxima'],
                'altura_minima': datos['altura_minima'],
                'num_islas': datos['num_islas'],
                'densidad_islas': datos['densidad_islas'],
                'cobertura_fraccional': datos['cobertura_fraccional'],
                'energia_difusion': datos['energia_difusion'],
                'tasa_eventos': datos['tasa_eventos'],
                'energia_promedio': datos['energia_promedio'],
                'etapa': datos['etapa']
            }
            
            # Añadir coberturas por umbral
            for umbral in [10, 20, 30, 60]:
                fila[f'cobertura_{umbral}'] = datos['coberturas_umbrales'].get(f'cobertura_{umbral}', 0)
            
            writer.writerow(fila)
    
    print(f"Resultados guardados en: {csv_path}")

def generar_reporte_texto(todos_datos):
    """Generar un reporte de texto con los resultados promedios"""
    if not todos_datos:
        print("No hay datos para generar reporte de texto")
        return
    
    # Separar datos por material
    datos_fe = [d for d in todos_datos if d['material'] == 'Fe']
    datos_cr = [d for d in todos_datos if d['material'] == 'Cr']
    
    # Calcular estadísticas para Fe
    estadisticas_fe = {}
    if datos_fe:
        estadisticas_fe = {
            'rugosidad': {'promedio': np.mean([d['rugosidad'] for d in datos_fe]), 
                         'desviacion': np.std([d['rugosidad'] for d in datos_fe])},
            'espesor_promedio': {'promedio': np.mean([d['espesor_promedio'] for d in datos_fe]), 
                                'desviacion': np.std([d['espesor_promedio'] for d in datos_fe])},
            'altura_maxima': {'promedio': np.mean([d['altura_maxima'] for d in datos_fe]), 
                             'desviacion': np.std([d['altura_maxima'] for d in datos_fe])},
            'etapas': {'nucleacion': sum(1 for d in datos_fe if d['etapa'] == 'nucleacion') / len(datos_fe) * 100,
                      'coalescencia': sum(1 for d in datos_fe if d['etapa'] == 'coalescencia') / len(datos_fe) * 100,
                      'crecimiento': sum(1 for d in datos_fe if d['etapa'] == 'crecimiento') / len(datos_fe) * 100}
        }
        
        # Calcular promedios para coberturas por umbral
        for umbral in [10, 20, 30, 60]:
            key = f'cobertura_{umbral}'
            valores = [d['coberturas_umbrales'].get(key, 0) for d in datos_fe]
            estadisticas_fe[key] = {
                'promedio': np.mean(valores),
                'desviacion': np.std(valores)
            }
    
    # Calcular estadísticas para Cr
    estadisticas_cr = {}
    if datos_cr:
        estadisticas_cr = {
            'rugosidad': {'promedio': np.mean([d['rugosidad'] for d in datos_cr]), 
                         'desviacion': np.std([d['rugosidad'] for d in datos_cr])},
            'espesor_promedio': {'promedio': np.mean([d['espesor_promedio'] for d in datos_cr]), 
                                'desviacion': np.std([d['espesor_promedio'] for d in datos_cr])},
            'altura_maxima': {'promedio': np.mean([d['altura_maxima'] for d in datos_cr]), 
                             'desviacion': np.std([d['altura_maxima'] for d in datos_cr])},
            'etapas': {'nucleacion': sum(1 for d in datos_cr if d['etapa'] == 'nucleacion') / len(datos_cr) * 100,
                      'coalescencia': sum(1 for d in datos_cr if d['etapa'] == 'coalescencia') / len(datos_cr) * 100,
                      'crecimiento': sum(1 for d in datos_cr if d['etapa'] == 'crecimiento') / len(datos_cr) * 100}
        }
        
        # Calcular promedios para coberturas por umbral
        for umbral in [10, 20, 30, 60]:
            key = f'cobertura_{umbral}'
            valores = [d['coberturas_umbrales'].get(key, 0) for d in datos_cr]
            estadisticas_cr[key] = {
                'promedio': np.mean(valores),
                'desviacion': np.std(valores)
            }
    
    # Generar reporte de texto
    reporte_path = os.path.join(output_dir, 'reporte_resultados.txt')
    
    with open(reporte_path, 'w') as f:
        f.write("RESULTADOS DEL ANÁLISIS DE SIMULACIONES\n")
        f.write("========================================\n\n")
        
        if datos_fe:
            f.write("Resultados para Fe:\n")
            f.write(f"Promedio - Rugosidad: {estadisticas_fe['rugosidad']['promedio']:.4f} ± {estadisticas_fe['rugosidad']['desviacion']:.4f}\n")
            f.write(f"Promedio - Espesor: {estadisticas_fe['espesor_promedio']['promedio']:.4f} ± {estadisticas_fe['espesor_promedio']['desviacion']:.4f}\n")
            f.write(f"Promedio - Altura máxima: {estadisticas_fe['altura_maxima']['promedio']:.4f} ± {estadisticas_fe['altura_maxima']['desviacion']:.4f}\n")
            
            for umbral in [10, 20, 30, 60]:
                key = f'cobertura_{umbral}'
                f.write(f"Promedio - Cobertura (> {umbral}): {estadisticas_fe[key]['promedio']:.2f}% ± {estadisticas_fe[key]['desviacion']:.2f}%\n")
            
            f.write("Distribución de etapas (%):\n")
            f.write(f"  nucleacion: {estadisticas_fe['etapas']['nucleacion']:.2f}%\n")
            f.write(f"  coalescencia: {estadisticas_fe['etapas']['coalescencia']:.2f}%\n")
            f.write(f"  crecimiento: {estadisticas_fe['etapas']['crecimiento']:.2f}%\n\n")
        
        if datos_cr:
            f.write("Resultados para Cr:\n")
            f.write(f"Promedio - Rugosidad: {estadisticas_cr['rugosidad']['promedio']:.4f} ± {estadisticas_cr['rugosidad']['desviacion']:.4f}\n")
            f.write(f"Promedio - Espesor: {estadisticas_cr['espesor_promedio']['promedio']:.4f} ± {estadisticas_cr['espesor_promedio']['desviacion']:.4f}\n")
            f.write(f"Promedio - Altura máxima: {estadisticas_cr['altura_maxima']['promedio']:.4f} ± {estadisticas_cr['altura_maxima']['desviacion']:.4f}\n")
            
            for umbral in [10, 20, 30, 60]:
                key = f'cobertura_{umbral}'
                f.write(f"Promedio - Cobertura (> {umbral}): {estadisticas_cr[key]['promedio']:.2f}% ± {estadisticas_cr[key]['desviacion']:.2f}%\n")
            
            f.write("Distribución de etapas (%):\n")
            f.write(f"  nucleacion: {estadisticas_cr['etapas']['nucleacion']:.2f}%\n")
            f.write(f"  coalescencia: {estadisticas_cr['etapas']['coalescencia']:.2f}%\n")
            f.write(f"  crecimiento: {estadisticas_cr['etapas']['crecimiento']:.2f}%\n\n")
        
        f.write(f"Total de simulaciones Fe: {len(datos_fe)}\n")
        f.write(f"Total de simulaciones Cr: {len(datos_cr)}\n")
        f.write(f"Total de simulaciones: {len(todos_datos)}\n")
    
    print(f"Reporte de texto generado en: {reporte_path}")

def main():
    """Función principal"""
    print("Iniciando análisis de simulaciones...")
    
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
            create_visualizations(archivo, resultados)
    
    if not todos_datos:
        print("No se pudieron procesar los archivos. Verifica el formato de los archivos .dat")
        return
    
    # Generar gráficos comparativos y promedios
    generar_graficos_comparativos(todos_datos)
    generar_graficos_promedio(todos_datos)
    
    # Guardar resultados en CSV
    guardar_resultados_csv(todos_datos)
    
    # Generar reporte de texto
    generar_reporte_texto(todos_datos)
    
    print(f"Proceso completado. Todos los resultados se han guardado en la carpeta '{output_dir}'")

if __name__ == "__main__":
    main()
