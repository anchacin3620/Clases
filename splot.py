import numpy as np
import matplotlib.pyplot as plt
import os
import re
from mpl_toolkits.mplot3d import Axes3D
from scipy import ndimage

# Definir umbrales para identificar etapas - AJUSTADOS para capturar etapas iniciales
UMBRALES_ETAPAS = {
    'nucleacion': 1.5,        # Reducido para capturar mejor la nucleación
    'coalescencia': (1.5, 5.0), # Ajustado para capturar mejor la coalescencia
    'crecimiento': 5.0        # Reducido para reflejar el menor número de eventos
}

def analizar_datos(archivo):
    """Analiza un archivo .dat generado por el código unificado y genera visualizaciones."""
    try:
        # Leer todo el contenido del archivo
        with open(archivo, 'r') as f:
            lines = f.readlines()
        
        # Separar encabezado y datos
        encabezado = []
        datos_lineas = []
        
        for line in lines:
            if line.strip().startswith('#'):
                encabezado.append(line.strip())
            else:
                datos_lineas.append(line.strip())
        
        # Procesar datos numéricos
        if not datos_lineas:
            print(f"Error: {archivo} no contiene datos válidos después del encabezado.")
            return None
        
        # Convertir datos a array numpy
        data = []
        for line in datos_lineas:
            if line:  # Asegurarse de que la línea no esté vacía
                partes = line.split()
                if len(partes) >= 3:  # Debe tener al menos 3 columnas
                    data.append([float(partes[0]), float(partes[1]), float(partes[2])])
        
        data = np.array(data)
        
        if data.size == 0:
            print(f"Error: {archivo} no contiene datos válidos después del encabezado.")
            return None
        
        # Determinar tamaño de la rejilla
        coordenadas_x = np.unique(data[:, 0])
        coordenadas_y = np.unique(data[:, 1])
        n_x = len(coordenadas_x)
        n_y = len(coordenadas_y)
        
        # Crear matriz de altura
        height = np.zeros((n_x, n_y))
        for punto in data:
            i = int(punto[0]) - 1  # Convertir a índice 0-based
            j = int(punto[1]) - 1  # Convertir a índice 0-based
            if 0 <= i < n_x and 0 <= j < n_y:
                height[i, j] = punto[2]
        
        # Calcular estadísticas básicas
        rugosidad = np.std(height)
        espesor_promedio = np.mean(height)
        max_altura = np.max(height)
        min_altura = np.min(height)
        
        # Extraer métricas del encabezado
        params = {}
        for line in encabezado:
            if '# Densidad de islas:' in line:
                params['densidad_islas'] = float(line.split(':')[1].strip())
            elif '# Tamaño promedio de islas:' in line:
                params['tamano_promedio_islas'] = float(line.split(':')[1].strip())
            elif '# Cobertura fraccional:' in line:
                params['cobertura_fraccional'] = float(line.split(':')[1].strip())
            elif '# Anisotropía superficial:' in line:
                params['anisotropia_superficial'] = float(line.split(':')[1].strip())
            elif '# Longitud de difusión promedio:' in line:
                params['longitud_difusion_promedio'] = float(line.split(':')[1].strip())
            elif '# Número de islas:' in line:
                params['num_islas'] = int(line.split(':')[1].strip())
            elif '# Átomos en superficie:' in line:
                params['atomos_superficie'] = int(line.split(':')[1].strip())
            elif '# Energía de difusión (eV):' in line:
                params['energia_difusion'] = float(line.split(':')[1].strip())
            elif '# Tasa de eventos:' in line:
                params['tasa_eventos'] = float(line.split(':')[1].strip())
            elif '# Energía promedio (eV):' in line:
                params['energia_promedio'] = float(line.split(':')[1].strip())
            elif '# Número de eventos:' in line:
                params['num_eventos'] = int(line.split(':')[1].strip())
            elif '# Rugosidad RMS (capas):' in line:
                params['rugosidad_rms'] = float(line.split(':')[1].strip())
            elif '# Tamaño de la rejilla:' in line:
                params['tamano_rejilla'] = int(line.split(':')[1].strip())
        
        # Si no se encuentran parámetros en el encabezado, intentar extraerlos del nombre del archivo
        nombre_archivo = os.path.basename(archivo)
        if 'energia_difusion' not in params:
            match_ed = re.search(r'ed([\d.]+)', nombre_archivo)
            if match_ed:
                params['energia_difusion'] = float(match_ed.group(1))
        
        if 'tasa_eventos' not in params:
            match_te = re.search(r'te([\d.]+)', nombre_archivo)
            if match_te:
                params['tasa_eventos'] = float(match_te.group(1))
        
        if 'energia_promedio' not in params:
            match_ep = re.search(r'ep([\d.]+)', nombre_archivo)
            if match_ep:
                params['energia_promedio'] = float(match_ep.group(1))
        
        match_s = re.search(r'_s(\d+)', nombre_archivo)
        if match_s:
            params['snapshot'] = int(match_s.group(1))
        
        match_c = re.search(r'_c(\d+)', nombre_archivo)
        if match_c:
            params['corrida'] = int(match_c.group(1))
        
        # Calcular métricas adicionales si no están en el encabezado
        if 'num_islas' not in params:
            # Calcular densidad de islas (número de islas por unidad de área)
            estructura = ndimage.generate_binary_structure(2, 2)
            mascara = height > 2  # Solo material depositado, excluyendo sustrato
            islas, num_islas = ndimage.label(mascara, structure=estructura)
            params['num_islas'] = num_islas
            params['densidad_islas'] = num_islas / (n_x * n_y)
        
        if 'cobertura_fraccional' not in params:
            # Calcular cobertura fraccional
            params['cobertura_fraccional'] = np.sum(height > 2) / (n_x * n_y) * 100
        
        # Calcular cobertura por umbrales
        umbrales = [2, 3, 5, 10]  # Umbrales ajustados para material depositado
        coberturas = {}
        for umbral in umbrales:
            areas_depositadas = height > umbral
            cobertura = np.sum(areas_depositadas) / (n_x * n_y) * 100
            coberturas[f'cobertura_{umbral}'] = cobertura
        
        # Preparar resultados para consola
        resultados_consola = [f"\nAnálisis de {archivo}:",
                             f"Rugosidad (desviación estándar): {rugosidad:.4f}",
                             f"Espesor promedio: {espesor_promedio:.4f}",
                             f"Altura máxima: {max_altura:.4f}",
                             f"Altura mínima: {min_altura:.4f}",
                             f"Número de islas: {params.get('num_islas', 'N/A')}",
                             f"Densidad de islas: {params.get('densidad_islas', 'N/A'):.6f}",
                             f"Tamaño promedio de islas: {params.get('tamano_promedio_islas', 'N/A'):.2f}",
                             f"Cobertura fraccional: {params.get('cobertura_fraccional', 'N/A'):.2f}%",
                             f"Anisotropía superficial: {params.get('anisotropia_superficial', 'N/A'):.4f}",
                             f"Longitud de difusión promedio: {params.get('longitud_difusion_promedio', 'N/A'):.4f}"]
        
        for umbral in umbrales:
            # CORRECCIÓN: Se eliminó el corchete extra en la f-string
            resultados_consola.append(f"Cobertura de islas (> {umbral}): {coberturas[f'cobertura_{umbral}']:.2f}%")
        
        # Visualización 2D
        plt.figure(figsize=(10, 8))
        plt.imshow(height, cmap='viridis', interpolation='nearest', vmin=0, vmax=max_altura)
        plt.colorbar(label='Altura (capas)')
        plt.title(f'Mapa 2D de Altura - {os.path.basename(archivo)}')
        plt.xlabel('Coordenada x (sitios)')
        plt.ylabel('Coordenada y (sitios)')
        
        carpeta_graficos = './graficos'
        if not os.path.exists(carpeta_graficos):
            os.makedirs(carpeta_graficos)
        nombre_base = os.path.basename(archivo).replace(".dat", "")
        plt.savefig(os.path.join(carpeta_graficos, f'mapa_2d_{nombre_base}.png'))
        plt.close()
        
        # Visualización 3D (solo para snapshots con menos de 10 capas promedio para mejor visualización)
        if espesor_promedio < 10:
            x = np.arange(n_x)
            y = np.arange(n_y)
            X, Y = np.meshgrid(x, y)
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X, Y, height, cmap='viridis', linewidth=0, antialiased=False, vmin=0, vmax=max_altura)
            ax.set_zlim(0, max_altura)
            ax.set_xlabel('Coordenada x (sitios)')
            ax.set_ylabel('Coordenada y (sitios)')
            ax.set_zlabel('Altura (capas)')
            ax.view_init(30, 45)
            plt.title(f'Mapa 3D de Altura - {os.path.basename(archivo)}')
            plt.colorbar(surf, shrink=0.5, aspect=5)
            plt.savefig(os.path.join(carpeta_graficos, f'mapa_3d_{nombre_base}.png'))
            plt.close()
        
        return {
            'height': height,
            'rugosidad': rugosidad, 
            'espesor_promedio': espesor_promedio,
            'altura_maxima': max_altura, 
            'altura_minima': min_altura,
            'n_x': n_x,
            'n_y': n_y,
            **params, 
            **coberturas, 
            'resultados_consola': resultados_consola
        }
    except Exception as e:
        print(f"Error al procesar {archivo}: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

def identificar_etapa(datos):
    """Identifica la etapa de crecimiento considerando espesor, cobertura y densidad de islas."""
    espesor = datos['espesor_promedio']
    cobertura = datos.get('cobertura_fraccional', 0)  # Usar cobertura fraccional del encabezado si está disponible
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

def analizar_corridas(directorio):
    """Analiza múltiples corridas and calcula estadísticas."""
    datos_totales = {'fe': [], 'cr': []}
    params = {'energia_difusion': [], 'tasa_eventos': [], 'energia_promedio': [], 'num_eventos': []}
    
    print(f"Buscando archivos en {directorio}...")
    for archivo in os.listdir(directorio):
        if archivo.endswith('.dat'):
            archivo_path = os.path.join(directorio, archivo)
            # Determinar material basado en el nombre del archivo
            if 'fe_' in archivo.lower():
                material = 'fe'
            elif 'cr_' in archivo.lower():
                material = 'cr'
            else:
                continue  # Saltar archivos que no son de Fe o Cr
                
            print(f"Procesando archivo: {archivo}")
            datos = analizar_datos(archivo_path)
            if datos is not None:
                datos['etapa'] = identificar_etapa(datos)
                datos_totales[material].append(datos)
                for key in ['energia_difusion', 'tasa_eventos', 'energia_promedio', 'num_eventos']:
                    if key in datos and datos[key] is not None:
                        params[key].append(datos[key])
            else:
                print(f"Omitiendo {archivo} debido a error de procesamiento.")
    
    resultados = {}
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            # Identificar todas las claves numéricas
            claves_numericas = set()
            for d in datos_totales[mat]:
                for k, v in d.items():
                    if k not in ['etapa', 'resultados_consola', 'height'] and isinstance(v, (int, float)):
                        claves_numericas.add(k)
            
            resultados[mat] = {
                'promedio': {k: np.mean([d[k] for d in datos_totales[mat] if d.get(k) is not None]) 
                            for k in claves_numericas},
                'desviacion': {k: np.std([d[k] for d in datos_totales[mat] if d.get(k) is not None]) 
                              for k in claves_numericas},
                'etapas_distribucion': {etapa: sum(1 for d in datos_totales[mat] if d['etapa'] == etapa) 
                                      / len(datos_totales[mat]) * 100 for etapa in ['nucleacion', 'coalescencia', 'crecimiento']}
            }
        else:
            resultados[mat] = {'promedio': {}, 'desviacion': {}, 'etapas_distribucion': {}}
            print(f"No se encontraron datos válidos para {mat.upper()}.")
    
    return datos_totales, resultados, params

def graficar_resultados(datos_totales, resultados, params):
    """Genera gráficos para comparar resultados, incluyendo histogramas adaptativos."""
    plt.figure(figsize=(20, 15))
    
    # Obtener el número real de snapshots para cada material
    n_fe = len(datos_totales['fe']) if datos_totales['fe'] else 0
    n_cr = len(datos_totales['cr']) if datos_totales['cr'] else 0
    
    # 1. Histograma Comparativo de Rugosidad (Fe vs. Cr)
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 1)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
                n_bins = min(15, max(5, len(rugosidades)//2))
                plt.hist(rugosidades, bins=n_bins, alpha=0.5, label=f'{mat.upper()} ({len(rugosidades)} snapshots)', density=True)
                plt.axvline(np.mean(rugosidades), color='k', linestyle='--', linewidth=1, label=f'Promedio {mat.upper()}')
        plt.xlabel('Valor de Rugosidad (capas)')
        plt.ylabel('Densidad de Probabilidad')
        plt.title(f'Distribución de la Rugosidad - {n_fe} snapshots Fe, {n_cr} snapshots Cr')
        plt.legend()
    
    # 2. Evolución de la rugosidad con el número de eventos
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 2)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                # Ordenar por número de eventos
                datos_ordenados = sorted(datos_totales[mat], key=lambda x: x.get('num_eventos', 0))
                eventos = [d.get('num_eventos', 0) for d in datos_ordenados]
                rugosidades = [d['rugosidad'] for d in datos_ordenados]
                plt.plot(eventos, rugosidades, 'o-', label=mat.upper(), alpha=0.7)
        plt.xlabel('Número de Eventos')
        plt.ylabel('Rugosidad (capas)')
        plt.title('Evolución de la Rugosidad con Eventos')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 3. Rugosidad vs Energía de difusión
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 3)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                energias = [d['energia_difusion'] for d in datos_totales[mat] if d.get('energia_difusion') is not None]
                rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
                if len(energias) > 0:
                    plt.scatter(energias, rugosidades, label=mat.upper(), alpha=0.7, s=80)
                    
                    if len(energias) > 1:
                        z = np.polyfit(energias, rugosidades, 1)
                        p = np.poly1d(z)
                        plt.plot(energias, p(energias), linestyle='--', alpha=0.5, linewidth=1)
                        
        plt.xlabel('Energía de difusión (eV)')
        plt.ylabel('Rugosidad (capas)')
        plt.title('Rugosidad vs Energía de difusión')
        plt.legend()
    
    # 4. Rugosidad vs Tasa de Eventos
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 4)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                tasas = [d['tasa_eventos'] for d in datos_totales[mat] if d.get('tasa_eventos') is not None]
                rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
                if len(tasas) > 0:
                    plt.scatter(tasas, rugosidades, label=mat.upper(), alpha=0.7, s=80)
                    
                    if len(tasas) > 1:
                        z = np.polyfit(tasas, rugosidades, 1)
                        p = np.poly1d(z)
                        plt.plot(tasas, p(tasas), linestyle='--', alpha=0.5, linewidth=1)
                        
        plt.xlabel('Tasa de Eventos')
        plt.ylabel('Rugosidad (capas)')
        plt.title('Rugosidad vs Tasa de Eventos')
        plt.legend()
    
    # 5. Rugosidad vs Energía Promedio
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 5)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                energias_prom = [d['energia_promedio'] for d in datos_totales[mat] if d.get('energia_promedio') is not None]
                rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
                if len(energias_prom) > 0:
                    plt.scatter(energias_prom, rugosidades, label=mat.upper(), alpha=0.7, s=80)
                    
                    if len(energias_prom) > 1:
                        z = np.polyfit(energias_prom, rugosidades, 1)
                        p = np.poly1d(z)
                        plt.plot(energias_prom, p(energias_prom), linestyle='--', alpha=0.5, linewidth=1)
                        
        plt.xlabel('Energía Promedio (eV)')
        plt.ylabel('Rugosidad (capas)')
        plt.title('Rugosidad vs Energía Promedio')
        plt.legend()
    
    # 6. Evolución del espesor con eventos
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 6)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                datos_ordenados = sorted(datos_totales[mat], key=lambda x: x.get('num_eventos', 0))
                eventos = [d.get('num_eventos', 0) for d in datos_ordenados]
                espesores = [d['espesor_promedio'] for d in datos_ordenados]
                plt.plot(eventos, espesores, 'o-', label=mat.upper(), alpha=0.7)
        plt.xlabel('Número de Eventos')
        plt.ylabel('Espesor Promedio (capas)')
        plt.title('Evolución del Espesor con Eventos')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 7. Densidad de islas vs espesor
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 7)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                espesores = [d['espesor_promedio'] for d in datos_totales[mat]]
                densidades = [d.get('densidad_islas', 0) for d in datos_totales[mat]]
                plt.scatter(espesores, densidades, label=mat.upper(), alpha=0.7, s=80)
        plt.xlabel('Espesor Promedio (capas)')
        plt.ylabel('Densidad de Islas')
        plt.title('Densidad de Islas vs Espesor')
        plt.legend()
        plt.yscale('log')  # Escala logarítmica para mejor visualización
    
    # 8. Cobertura fraccional vs espesor
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 8)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                espesores = [d['espesor_promedio'] for d in datos_totales[mat]]
                coberturas = [d.get('cobertura_fraccional', 0) for d in datos_totales[mat]]
                plt.scatter(espesores, coberturas, label=mat.upper(), alpha=0.7, s=80)
        plt.xlabel('Espesor Promedio (capas)')
        plt.ylabel('Cobertura Fraccional (%)')
        plt.title('Cobertura Fraccional vs Espesor')
        plt.legend()
    
    # 9. Tamaño promedio de islas vs espesor
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 9)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                espesores = [d['espesor_promedio'] for d in datos_totales[mat]]
                tamanos = [d.get('tamano_promedio_islas', 0) for d in datos_totales[mat]]
                plt.scatter(espesores, tamanos, label=mat.upper(), alpha=0.7, s=80)
        plt.xlabel('Espesor Promedio (capas)')
        plt.ylabel('Tamaño Promedio de Islas')
        plt.title('Tamaño Promedio de Islas vs Espesor')
        plt.legend()
    
    # 10. Anisotropía superficial vs espesor
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 10)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                espesores = [d['espesor_promedio'] for d in datos_totales[mat]]
                anisotropias = [d.get('anisotropia_superficial', 0) for d in datos_totales[mat]]
                plt.scatter(espesores, anisotropias, label=mat.upper(), alpha=0.7, s=80)
        plt.xlabel('Espesor Promedio (capas)')
        plt.ylabel('Anisotropía Superficial')
        plt.title('Anisotropía Superficial vs Espesor')
        plt.legend()
    
    # 11. Longitud de difusión vs energía de difusión
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 11)
        for mat in ['fe', 'cr']:
            if datos_totales[mat]:
                energias = [d['energia_difusion'] for d in datos_totales[mat] if d.get('energia_difusion') is not None]
                longitudes = [d.get('longitud_difusion_promedio', 0) for d in datos_totales[mat]]
                if len(energias) > 0:
                    plt.scatter(energias, longitudes, label=mat.upper(), alpha=0.7, s=80)
        plt.xlabel('Energía de difusión (eV)')
        plt.ylabel('Longitud de Difusión Promedio')
        plt.title('Longitud de Difusión vs Energía de Difusión')
        plt.legend()
    
    # 12. Etapas de crecimiento
    if n_fe > 0 or n_cr > 0:
        plt.subplot(3, 4, 12)
        for mat in ['fe', 'cr']:
            if resultados[mat]['etapas_distribucion']:
                etapas = list(resultados[mat]['etapas_distribucion'].keys())
                porcentajes = [resultados[mat]['etapas_distribucion'][etapa] for etapa in etapas]
                n_snapshots = len(datos_totales[mat])
                plt.bar([f'{mat.upper()}_{etapa}' for etapa in etapas],
                       porcentajes, label=f'{mat.upper()} ({n_snapshots} snapshots)', alpha=0.7)
        
        plt.xlabel('Etapa')
        plt.ylabel('Porcentaje (%)')
        plt.title('Distribución de etapas de crecimiento')
        plt.legend()
        plt.xticks(rotation=45)
    
    # Ajustar diseño y guardar
    plt.tight_layout()
    
    # Guardar gráfico comparativo
    carpeta_graficos = './graficos'
    if not os.path.exists(carpeta_graficos):
        os.makedirs(carpeta_graficos)
    plt.savefig(os.path.join(carpeta_graficos, 'grafico_comparativo_completo.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Función principal que ejecuta el análisis completo."""
    # CORRECCIÓN: Cambiar el directorio a './simulaciones' en lugar del directorio actual
    directorio = './simulaciones'
    if not os.path.exists(directorio):
        print(f"Directorio {directorio} no encontrado. Creándolo...")
        os.makedirs(directorio)
        print(f"Por favor, coloca los archivos .dat generados por Fortran en la carpeta {directorio} y ejecuta nuevamente.")
        return
    
    datos_totales, resultados, params = analizar_corridas(directorio)
    
    # Almacena resultados en archivo .txt
    archivo_resultados = 'resultados_analisis.txt'
    with open(archivo_resultados, 'w') as f:
        for mat in ['fe', 'cr']:
            if resultados[mat]['promedio']:
                f.write(f"\nResultados para {mat.upper()}:\n")
                f.write(f"Promedio - Rugosidad: {resultados[mat]['promedio'].get('rugosidad', 'N/A'):.4f} ± "
                        f"{resultados[mat]['desviacion'].get('rugosidad', 'N/A'):.4f}\n")
                f.write(f"Promedio - Espesor: {resultados[mat]['promedio'].get('espesor_promedio', 'N/A'):.4f} ± "
                        f"{resultados[mat]['desviacion'].get('espesor_promedio', 'N/A'):.4f}\n")
                f.write(f"Promedio - Altura máxima: {resultados[mat]['promedio'].get('altura_maxima', 'N/A'):.4f} ± "
                        f"{resultados[mat]['desviacion'].get('altura_maxima', 'N/A'):.4f}\n")
                f.write(f"Promedio - Número de islas: {resultados[mat]['promedio'].get('num_islas', 'N/A'):.2f} ± "
                        f"{resultados[mat]['desviacion'].get('num_islas', 'N/A'):.2f}\n")
                f.write(f"Promedio - Densidad de islas: {resultados[mat]['promedio'].get('densidad_islas', 'N/A'):.6f} ± "
                        f"{resultados[mat]['desviacion'].get('densidad_islas', 'N/A'):.6f}\n")
                f.write(f"Promedio - Tamaño promedio de islas: {resultados[mat]['promedio'].get('tamano_promedio_islas', 'N/A'):.2f} ± "
                        f"{resultados[mat]['desviacion'.get('tamano_promedio_islas', 'N/A'):.2f}\n")
                f.write(f"Promedio - Cobertura fraccional: {resultados[mat]['promedio'].get('cobertura_fraccional', 'N/A'):.2f}% ± "
                        f"{resultados[mat]['desviacion'].get('cobertura_fraccional', 'N/A'):.2f}%\n")
                f.write(f"Promedio - Anisotropía superficial: {resultados[mat]['promedio'].get('anisotropia_superficial', 'N/A'):.4f} ± "
                        f"{resultados[mat]['desviacion'].get('anisotropia_superficial', 'N/A'):.4f}\n")
                f.write(f"Promedio - Longitud de difusión promedio: {resultados[mat]['promedio'].get('longitud_difusion_promedio', 'N/A'):.4f} ± "
                        f"{resultados[mat]['desviacion'].get('longitud_difusion_promedio', 'N/A'):.4f}\n")
                
                for umbral in [2, 3, 5, 10]:
                    f.write(f"Promedio - Cobertura (> {umbral}): "
                            f"{resultados[mat]['promedio'].get(f'cobertura_{umbral}', 'N/A'):.2f}% ± "
                            f"{resultados[mat]['desviacion'].get(f'cobertura_{umbral}', 'N/A'):.2f}%\n")
                
                f.write("Distribución de etapas (%):\n")
                for etapa, porcentaje in resultados[mat]['etapas_distribucion'].items():
                    f.write(f"  {etapa}: {porcentaje:.2f}%\n")
            else:
                f.write(f"\nNo se encontraron datos válidos para {mat.upper()}.\n")
        
        if not any(resultados[mat]['promedio'] for mat in ['fe', 'cr']):
            f.write("No se encontraron datos para graficar.\n")
    
    # Imprime resultados en consola
    print(f"Resultados guardados en {archivo_resultados}")
    for mat in ['fe', 'cr']:
        if resultados[mat]['promedio']:
            print(f"\nResultados para {mat.upper()}:")
            print(f"Promedio - Rugosidad: {resultados[mat]['promedio'].get('rugosidad', 'N/A'):.4f} ± "
                  f"{resultados[mat]['desviacion'].get('rugosidad', 'N/A'):.4f}")
            print(f"Promedio - Espesor: {resultados[mat]['promedio'].get('espesor_promedio', 'N/A'):.4f} ± "
                  f"{resultados[mat]['desviacion'].get('espesor_promedio', 'N/A'):.4f}")
            print(f"Promedio - Altura máxima: {resultados[mat]['promedio'].get('altura_maxima', 'N/A'):.4f} ± "
                  f"{resultados[mat]['desviacion'].get('altura_maxima', 'N/A'):.4f}")
            print(f"Promedio - Número de islas: {resultados[mat]['promedio'].get('num_islas', 'N/A'):.2f} ± "
                  f"{resultados[mat]['desviacion'].get('num_islas', 'N/A'):.2f}")
            print(f"Promedio - Densidad de islas: {resultados[mat]['promedio'].get('densidad_islas', 'N/A'):.6f} ± "
                  f"{resultados[mat]['desviacion'].get('densidad_islas', 'N/A'):.6f}")
            print(f"Promedio - Tamaño promedio de islas: {resultados[mat]['promedio'].get('tamano_promedio_islas', 'N/A'):.2f} ± "
                  f"{resultados[mat]['desviacion'].get('tamano_promedio_islas', 'N/A'):.2f}")
            print(f"Promedio - Cobertura fraccional: {resultados[mat]['promedio'].get('cobertura_fraccional', 'N/A'):.2f}% ± "
                  f"{resultados[mat]['desviacion'].get('cobertura_fraccional', 'N/A'):.2f}%")
            print(f"Promedio - Anisotropía superficial: {resultados[mat]['promedio'].get('anisotropia_superficial', 'N/A'):.4f} ± "
                  f"{resultados[mat]['desviacion'].get('anisotropia_superficial', 'N/A'):.4f}")
            print(f"Promedio - Longitud de difusión promedio: {resultados[mat]['promedio'].get('longitud_difusion_promedio', 'N/A'):.4f} ± "
                  f"{resultados[mat]['desviacion'].get('longitud_difusion_promedio', 'N/A'):.4f}")
            
            for umbral in [2, 3, 5, 10]:
                print(f"Promedio - Cobertura (> {umbral}): "
                      f"{resultados[mat]['promedio'].get(f'cobertura_{umbral}', 'N/A'):.2f}% ± "
                      f"{resultados[mat]['desviacion'].get(f'cobertura_{umbral}', 'N/A'):.2f}%")
            
            print("Distribución de etapas (%):")
            for etapa, porcentaje in resultados[mat]['etapas_distribucion'].items():
                print(f"  {etapa}: {porcentaje:.2f}%")
        else:
            print(f"\nNo se encontraron datos válidos para {mat.upper()}.")
    
    # Genera gráficos si hay datos
    if any(resultados[mat]['promedio'] for mat in ['fe', 'cr']):
        graficar_resultados(datos_totales, resultados, params)
        print("Gráficos generados correctamente.")
    else:
        print("No se encontraron datos para graficar.")

if __name__ == "__main__":
    main()
