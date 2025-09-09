import numpy as np
import matplotlib.pyplot as plt
import os
import re
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from scipy.stats import gaussian_kde

# Definir umbrales para identificar etapas
UMBRALES_ETAPAS = {
    'nucleacion': 5,
    'coalescencia': (5, 15),
    'crecimiento': 15
}

def analizar_datos(archivo):
    """Analiza un archivo .dat generado por los códigos Fortran y genera visualizaciones."""
    try:
        # Contar líneas de encabezado que comienzan con '#'
        with open(archivo, 'r') as f:
            num_encabezado = sum(1 for line in f if line.strip().startswith('#'))
            f.seek(0)  # Regresa al inicio del archivo para leer el encabezado
            encabezado = [line.strip() for line in f if line.strip().startswith('#')]

        # Leer datos numéricos, omitiendo el encabezado
        data = np.loadtxt(archivo, skiprows=num_encabezado)
        if data.size == 0:
            print(f"Error: {archivo} no contiene datos válidos después del encabezado.")
            return None
        n = int(np.sqrt(len(data)))  # Tamaño de la rejilla (e.g., 300)
        height = data[:, 2].reshape(n, n)  # Extrae altura y reestructura

        # Calcular estadísticas
        rugosidad = np.std(height)
        espesor_promedio = np.mean(height)
        max_altura = np.max(height)
        min_altura = np.min(height)

        resultados_consola = [f"\nAnálisis de {archivo}:",
                             f"Rugosidad (desviación estándar): {rugosidad:.4f}",
                             f"Espesor promedio: {espesor_promedio:.4f}",
                             f"Altura máxima: {max_altura:.4f}",
                             f"Altura mínima: {min_altura:.4f}"]

        # Calcular cobertura de islas
        umbrales = [10, 20, 30, 60]
        coberturas = {}
        for umbral in umbrales:
            islas = height > umbral
            cobertura = np.sum(islas) / (n * n) * 100
            coberturas[f'cobertura_{umbral}'] = cobertura
            resultados_consola.append(f"Cobertura de islas (> {umbral}): {cobertura:.2f}%")

        # Visualización 2D
        plt.figure(figsize=(10, 8))
        plt.imshow(height, cmap='viridis', interpolation='nearest', vmin=0, vmax=max_altura)
        plt.colorbar(label='Altura (capas)')
        plt.title(f'Mapa 2D de Altura - {archivo}')
        plt.xlabel('Coordenada x (sitios)')
        plt.ylabel('Coordenada y (sitios)')
        
        carpeta_graficos = './graficos'
        if not os.path.exists(carpeta_graficos):
            os.makedirs(carpeta_graficos)
        plt.savefig(os.path.join(carpeta_graficos, f'mapa_2d_{os.path.basename(archivo).replace(".dat", ".png")}'))
        plt.close()

        # Visualización 3D
        x = np.arange(n)
        y = np.arange(n)
        X, Y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X, Y, height, cmap='viridis', linewidth=0, antialiased=False, vmin=0, vmax=max_altura)
        ax.set_zlim(0, max_altura)
        ax.set_xlabel('Coordenada x (sitios)')
        ax.set_ylabel('Coordenada y (sitios)')
        ax.set_zlabel('Altura (capas)')
        ax.view_init(30, 45)
        plt.title(f'Mapa 3D de Altura - {archivo}')
        plt.colorbar(surf, shrink=0.5, aspect=5)
        plt.savefig(os.path.join(carpeta_graficos, f'mapa_3d_{os.path.basename(archivo).replace(".dat", ".png")}'))
        plt.close()

        # Extraer parámetros del encabezado con manejo de errores
        params = {}
        with open(archivo, 'r') as f:
            for line in f:
                if line.strip().startswith('# Energía de difusión'):
                    match = re.search(r'[-+]?\d*\.\d+|\d+', line)
                    params['energia_difusion'] = float(match.group()) if match else None
                elif line.strip().startswith('# Tasa de eventos'):
                    match = re.search(r'[-+]?\d*\.\d+|\d+', line)
                    params['tasa_eventos'] = float(match.group()) if match else None
                elif line.strip().startswith('# Energía promedio'):
                    match = re.search(r'[-+]?\d*\.\d+|\d+', line)
                    params['energia_promedio'] = float(match.group()) if match else None

        # Si no se encuentran parámetros en el encabezado, intentar extraerlos del nombre del archivo
        if not all(params.values()):
            nombre = os.path.basename(archivo).replace('.dat', '')
            match_ed = re.search(r'ed([\d.]+)', nombre)
            match_te = re.search(r'te([\d.]+)', nombre)
            match_ep = re.search(r'ep\s*([\d.]+)', nombre)
            if match_ed:
                params['energia_difusion'] = float(match_ed.group(1))
            if match_te:
                params['tasa_eventos'] = float(match_te.group(1))
            if match_ep:
                params['energia_promedio'] = float(match_ep.group(1))

        return {
            'height': height,
            'rugosidad': rugosidad, 'espesor_promedio': espesor_promedio,
            'altura_maxima': max_altura, 'altura_minima': min_altura,
            **coberturas, **{k: v for k, v in params.items() if v is not None}, 
            'resultados_consola': resultados_consola
        }
    except Exception as e:
        print(f"Error al procesar {archivo}: {str(e)}")
        return None

def identificar_etapa(espesor_promedio):
    """Identifica la etapa de crecimiento."""
    if espesor_promedio < UMBRALES_ETAPAS['nucleacion']:
        return 'nucleacion'
    elif UMBRALES_ETAPAS['coalescencia'][0] <= espesor_promedio <= UMBRALES_ETAPAS['coalescencia'][1]:
        return 'coalescencia'
    else:
        return 'crecimiento'

def analizar_corridas(directorio):
    """Analiza múltiples corridas y calcula estadísticas."""
    datos_totales = {'fe': [], 'cr': []}
    params = {'energia_difusion': [], 'tasa_eventos': [], 'energia_promedio': []}

    print(f"Buscando archivos en {directorio}...")
    for archivo in os.listdir(directorio):
        if archivo.endswith('.dat') and ('fe' in archivo.lower() or 'cr' in archivo.lower()):
            print(f"Detectado archivo: {archivo}")
            material = 'fe' if 'fe' in archivo.lower() else 'cr'
            datos = analizar_datos(os.path.join(directorio, archivo))
            if datos is not None:
                datos['etapa'] = identificar_etapa(datos['espesor_promedio'])
                datos_totales[material].append(datos)
                for key in ['energia_difusion', 'tasa_eventos', 'energia_promedio']:
                    if key in datos and datos[key] is not None:
                        params[key].append(datos[key])
            else:
                print(f"Omitiendo {archivo} debido a error de procesamiento.")

    resultados = {}
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            claves_numericas = [k for k in datos_totales[mat][0].keys() 
                              if k not in ['etapa', 'resultados_consola', 'height']]
            resultados[mat] = {
                'promedio': {k: np.mean([d[k] for d in datos_totales[mat] if d.get(k) is not None]) 
                            for k in claves_numericas},
                'desviacion': {k: np.std([d[k] for d in datos_totales[mat] if d.get(k) is not None]) 
                              for k in claves_numericas},
                'etapas_distribucion': {etapa: sum(1 for d in datos_totales[mat] if d['etapa'] == etapa) 
                                      / len(datos_totales[mat]) * 100 for etapa in UMBRALES_ETAPAS.keys()}
            }
        else:
            resultados[mat] = {'promedio': {}, 'desviacion': {}, 'etapas_distribucion': {}}
            print(f"No se encontraron datos válidos para {mat.upper()}.")

    return datos_totales, resultados, params

def graficar_resultados(datos_totales, resultados, params):
    """Genera gráficos para comparar resultados, incluyendo histogramas adaptativos."""
    plt.figure(figsize=(15, 10))
    
    # Obtener el número real de corridas para cada material
    n_fe = len(datos_totales['fe']) if datos_totales['fe'] else 0
    n_cr = len(datos_totales['cr']) if datos_totales['cr'] else 0
    
    # Histograma Comparativo de Rugosidad (Fe vs. Cr)
    plt.subplot(2, 2, 1)
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
            # Ajustar bins según el número de corridas
            n_bins = min(15, max(5, len(rugosidades)//2))
            plt.hist(rugosidades, bins=n_bins, alpha=0.5, label=f'{mat.upper()} ({len(rugosidades)} corridas)', density=True)
            plt.axvline(np.mean(rugosidades), color='k', linestyle='--', linewidth=1, label=f'Promedio {mat.upper()}')
    plt.xlabel('Valor de Rugosidad (capas)')
    plt.ylabel('Densidad de Probabilidad')
    plt.title(f'Distribución de la Rugosidad - {n_fe} corridas Fe, {n_cr} corridas Cr')
    plt.legend()

    # Histograma de Convergencia del Cromo (Cr) - Adaptativo
    plt.subplot(2, 2, 2)
    if datos_totales['cr']:
        n_corridas = len(datos_totales['cr'])
        
        # Definir puntos de muestra según el número de corridas
        if n_corridas >= 10:
            puntos_muestra = [2, 5, 10]
        elif n_corridas >= 5:
            puntos_muestra = [2, n_corridas]
        elif n_corridas >= 2:
            puntos_muestra = [2]
        else:
            puntos_muestra = [1]
        
        # Ajustar el máximo de puntos a mostrar
        puntos_muestra = puntos_muestra[:min(3, n_corridas)]
        
        for num_corridas in puntos_muestra:
            rugosidades = [d['rugosidad'] for d in datos_totales['cr'][:num_corridas]]
            n_bins = min(10, max(3, num_corridas))
            plt.hist(rugosidades, bins=n_bins, alpha=0.6, 
                    label=f'{num_corridas} Corrida{"s" if num_corridas > 1 else ""}', 
                    density=True, histtype='stepfilled')
            plt.axvline(np.mean(rugosidades), linestyle='--', linewidth=1.5, 
                       label=f'Promedio {num_corridas} corrida{"s" if num_corridas > 1 else ""}')
        
        plt.xlabel('Valor de Rugosidad (capas)')
        plt.ylabel('Densidad de Probabilidad')
        plt.title(f'Convergencia de la Rugosidad del Cromo ({n_corridas} corridas totales)')
        plt.legend()

    # Rugosidad vs Energía de difusión (por corrida)
    plt.subplot(2, 2, 3)
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            energias = [d['energia_difusion'] for d in datos_totales[mat] if d.get('energia_difusion') is not None]
            rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
            plt.scatter(energias, rugosidades, label=mat.upper(), alpha=0.7, s=80)
            
            # Añadir línea de tendencia si hay suficientes puntos
            if len(energias) > 1:
                z = np.polyfit(energias, rugosidades, 1)
                p = np.poly1d(z)
                plt.plot(energias, p(energias), linestyle='--', alpha=0.5, linewidth=1)
                
    plt.xlabel('Energía de difusión (eV)')
    plt.ylabel('Rugosidad')
    plt.title('Rugosidad vs Energía de difusión (por corrida)')
    plt.legend()

    # Etapas de crecimiento
    plt.subplot(2, 2, 4)
    for mat in ['fe', 'cr']:
        if resultados[mat]['etapas_distribucion']:
            etapas = list(UMBRALES_ETAPAS.keys())
            porcentajes = [resultados[mat]['etapas_distribucion'][etapa] for etapa in etapas]
            
            # Añadir etiqueta con número de corridas
            n_corridas_mat = len(datos_totales[mat])
            plt.bar([f'{mat.upper()}_{etapa}' for etapa in etapas],
                   porcentajes, label=f'{mat.upper()} ({n_corridas_mat} corridas)', alpha=0.7)
    
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
    plt.savefig(os.path.join(carpeta_graficos, 'grafico_comparativo_mejorado.png'), dpi=300, bbox_inches='tight')
    plt.close()

def graficar_histogramas_individuales(datos_totales):
    """Genera histogramas individuales para cada material y número de corridas."""
    for mat in ['fe', 'cr']:
        if not datos_totales[mat]:
            continue
            
        n_corridas = len(datos_totales[mat])
        rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
        
        plt.figure(figsize=(10, 6))
        n_bins = min(15, max(5, n_corridas))
        
        # Histograma con curva de densidad
        n, bins, patches = plt.hist(rugosidades, bins=n_bins, alpha=0.7, 
                                   density=True, label=f'Rugosidad {mat.upper()}')
        
        # Añadir línea de densidad
        if len(rugosidades) > 1:
            density = gaussian_kde(rugosidades)
            xs = np.linspace(min(rugosidades), max(rugosidades), 200)
            plt.plot(xs, density(xs), 'r-', linewidth=2, label='Densidad')
        
        # Líneas de media y mediana
        mean_rug = np.mean(rugosidades)
        median_rug = np.median(rugosidades)
        plt.axvline(mean_rug, color='green', linestyle='--', linewidth=2, 
                   label=f'Media: {mean_rug:.2f}')
        plt.axvline(median_rug, color='blue', linestyle='--', linewidth=2, 
                   label=f'Mediana: {median_rug:.2f}')
        
        plt.xlabel('Rugosidad (capas)')
        plt.ylabel('Densidad de Probabilidad')
        plt.title(f'Distribución de Rugosidad - {mat.upper()} ({n_corridas} corridas)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Guardar
        carpeta_graficos = './graficos'
        plt.savefig(os.path.join(carpeta_graficos, f'histograma_{mat}_{n_corridas}corridas.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()

def graficar_evolucion_rugosidad(datos_totales):
    """Genera gráfico de evolución de rugosidad por número de corridas."""
    plt.figure(figsize=(12, 6))
    
    for mat in ['fe', 'cr']:
        if not datos_totales[mat]:
            continue
            
        rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
        # Calcular media acumulativa
        medias_acumulativas = [np.mean(rugosidades[:i+1]) for i in range(len(rugosidades))]
        # Calcular desviación estándar acumulativa
        desviaciones_acumulativas = [np.std(rugosidades[:i+1]) for i in range(len(rugosidades))]
        
        corridas = range(1, len(rugosidades)+1)
        plt.errorbar(corridas, medias_acumulativas, yerr=desviaciones_acumulativas, 
                    fmt='-o', capsize=4, label=mat.upper(), alpha=0.8)
    
    plt.xlabel('Número de Corridas')
    plt.ylabel('Rugosidad Promedio (capas)')
    plt.title('Evolución de la Rugosidad con Aumento de Corridas')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Guardar
    carpeta_graficos = './graficos'
    plt.savefig(os.path.join(carpeta_graficos, 'evolucion_rugosidad.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()

def graficar_promedios(datos_totales, resultados):
    """Genera gráficos promedios de altura por material."""
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            # Calcular altura promedio
            heights = [d['height'] for d in datos_totales[mat]]
            altura_promedio = np.mean(heights, axis=0)
            n = altura_promedio.shape[0]

            # Visualización 2D promedio
            plt.figure(figsize=(10, 8))
            plt.imshow(altura_promedio, cmap='viridis', interpolation='nearest', vmin=0, vmax=np.max(altura_promedio))
            plt.colorbar(label='Altura promedio (capas)')
            plt.title(f'Mapa 2D de Altura Promedio - {mat.upper()}')
            plt.xlabel('Coordenada x (sitios)')
            plt.ylabel('Coordenada y (sitios)')
            plt.savefig(os.path.join('./graficos', f'mapa_2d_promedio_{mat}.png'))
            plt.close()

            # Visualización 3D promedio
            x = np.arange(n)
            y = np.arange(n)
            X, Y = np.meshgrid(x, y)
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X, Y, altura_promedio, cmap='viridis', linewidth=0, antialiased=False, vmin=0, vmax=np.max(altura_promedio))
            ax.set_zlim(0, np.max(altura_promedio))
            ax.set_xlabel('Coordenada x (sitios)')
            ax.set_ylabel('Coordenada y (sitios)')
            ax.set_zlabel('Altura promedio (capas)')
            ax.view_init(30, 45)
            plt.title(f'Mapa 3D de Altura Promedio - {mat.upper()}')
            plt.colorbar(surf, shrink=0.5, aspect=5)
            plt.savefig(os.path.join('./graficos', f'mapa_3d_promedio_{mat}.png'))
            plt.close()

def main():
    """Función principal que ejecuta el análisis completo."""
    directorio = './simulaciones'
    if not os.path.exists(directorio):
        print(f"Directorio {directorio} no encontrado. Creándolo...")
        os.makedirs(directorio)

    datos_totales, resultados, params = analizar_corridas(directorio)

    # Almacena resultados en archivo .txt
    archivo_resultados = 'resultados_analisis.txt'
    with open(archivo_resultados, 'w') as f:
        for mat in ['fe', 'cr']:
            if resultados[mat]['promedio']:
                f.write(f"\nResultados para {mat.upper()}:\n")
                f.write(f"Promedio - Rugosidad: {resultados[mat]['promedio']['rugosidad']:.4f} ± "
                        f"{resultados[mat]['desviacion']['rugosidad']:.4f}\n")
                f.write(f"Promedio - Espesor: {resultados[mat]['promedio']['espesor_promedio']:.4f} ± "
                        f"{resultados[mat]['desviacion']['espesor_promedio']:.4f}\n")
                f.write(f"Promedio - Altura máxima: {resultados[mat]['promedio']['altura_maxima']:.4f} ± "
                        f"{resultados[mat]['desviacion']['altura_maxima']:.4f}\n")
                for umbral in [10, 20, 30, 60]:
                    f.write(f"Promedio - Cobertura (> {umbral}): "
                            f"{resultados[mat]['promedio'][f'cobertura_{umbral}']:.2f}% ± "
                            f"{resultados[mat]['desviacion'][f'cobertura_{umbral}']:.2f}%\n")
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
            print(f"Promedio - Rugosidad: {resultados[mat]['promedio']['rugosidad']:.4f} ± "
                  f"{resultados[mat]['desviacion']['rugosidad']:.4f}")
            print(f"Promedio - Espesor: {resultados[mat]['promedio']['espesor_promedio']:.4f} ± "
                  f"{resultados[mat]['desviacion']['espesor_promedio']:.4f}")
            print(f"Promedio - Altura máxima: {resultados[mat]['promedio']['altura_maxima']:.4f} ± "
                  f"{resultados[mat]['desviacion']['altura_maxima']:.4f}")
            for umbral in [10, 20, 30, 60]:
                print(f"Promedio - Cobertura (> {umbral}): "
                      f"{resultados[mat]['promedio'][f'cobertura_{umbral}']:.2f}% ± "
                      f"{resultados[mat]['desviacion'][f'cobertura_{umbral}']:.2f}%")
            print("Distribución de etapas (%):")
            for etapa, porcentaje in resultados[mat]['etapas_distribucion'].items():
                print(f"  {etapa}: {porcentaje:.2f}%")
        else:
            print(f"\nNo se encontraron datos válidos para {mat.upper()}.")

    # Genera gráficos si hay datos
    if any(resultados[mat]['promedio'] for mat in ['fe', 'cr']):
        graficar_resultados(datos_totales, resultados, params)
        graficar_promedios(datos_totales, resultados)
        graficar_histogramas_individuales(datos_totales)
        graficar_evolucion_rugosidad(datos_totales)
        print("Gráficos generados correctamente.")
    else:
        print("No se encontraron datos para graficar.")

if __name__ == "__main__":
    main()
