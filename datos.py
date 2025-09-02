import numpy as np
import matplotlib.pyplot as plt
import os
import re
from mpl_toolkits.mplot3d import Axes3D

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
            print(f"Procesando {archivo} con {num_encabezado} líneas de encabezado.")
            f.seek(0)  # Regresa al inicio del archivo para leer el encabezado
            encabezado = [line.strip() for line in f if line.strip().startswith('#')]
            print("Contenido del encabezado:", encabezado)

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
            print(f"Parámetros extraídos del nombre {nombre}: {params}")

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
                print(f"Omitting {archivo} due to processing error.")

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
    """Genera gráficos para comparar resultados, incluyendo histogramas."""
    plt.figure(figsize=(15, 10))

    plt.subplot(2, 2, 1)
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
            plt.hist(rugosidades, bins=10, alpha=0.5, label=mat.upper(), density=True)
            plt.axvline(np.mean(rugosidades), color='k', linestyle='--', linewidth=1, 
                       label=f'Promedio {mat.upper()}')
    plt.xlabel('Valor de Rugosidad (capas)')
    plt.ylabel('Frecuencia de Ocurrencia')
    plt.title('Distribución de la Rugosidad (10 Corridas)')
    plt.legend()

    plt.subplot(2, 2, 2)
    if datos_totales['cr']:
        rugosidades_2 = [d['rugosidad'] for d in datos_totales['cr'][:2]]
        rugosidades_5 = [d['rugosidad'] for d in datos_totales['cr'][:5]]
        rugosidades_10 = [d['rugosidad'] for d in datos_totales['cr']]
        plt.hist(rugosidades_2, bins=10, alpha=0.5, label='2 Corridas', density=True)
        plt.hist(rugosidades_5, bins=10, alpha=0.5, label='5 Corridas', density=True)
        plt.hist(rugosidades_10, bins=10, alpha=0.5, label='10 Corridas', density=True)
        plt.axvline(np.mean(rugosidades_2), color='k', linestyle='--', linewidth=1, 
                   label='Promedio 2 Corridas')
        plt.axvline(np.mean(rugosidades_5), color='k', linestyle='--', linewidth=1, 
                   label='Promedio 5 Corridas')
        plt.axvline(np.mean(rugosidades_10), color='k', linestyle='--', linewidth=1, 
                   label='Promedio 10 Corridas')
    plt.xlabel('Valor de Rugosidad (capas)')
    plt.ylabel('Frecuencia de Ocurrencia')
    plt.title('Estabilización de la Distribución de la Rugosidad del Cromo')
    plt.legend()

    plt.subplot(2, 2, 3)
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            energias = [d['energia_difusion'] for d in datos_totales[mat] if d.get('energia_difusion') is not None]
            rugosidades = [d['rugosidad'] for d in datos_totales[mat]]
            plt.scatter(energias, rugosidades, label=mat.upper(), alpha=0.5)
    plt.xlabel('Energía de difusión (eV)')
    plt.ylabel('Rugosidad')
    plt.title('Rugosidad vs Energía de difusión (por corrida)')
    plt.legend()

    plt.subplot(2, 2, 4)
    for mat in ['fe', 'cr']:
        if resultados[mat]['etapas_distribucion']:
            plt.bar([f'{mat.upper()}_{etapa}' for etapa in UMBRALES_ETAPAS.keys()],
                   [resultados[mat]['etapas_distribucion'][etapa] for etapa in UMBRALES_ETAPAS.keys()],
                   label=mat.upper())
    plt.xlabel('Etapa')
    plt.ylabel('Porcentaje (%)')
    plt.title('Distribución de etapas de crecimiento')
    plt.legend()
    
    carpeta_graficos = './graficos'
    if not os.path.exists(carpeta_graficos):
        os.makedirs(carpeta_graficos)
    plt.savefig(os.path.join(carpeta_graficos, 'grafico_comparativo.png'))
    plt.close()

def graficar_promedios(datos_totales, resultados):
    """Genera gráficos promedios de altura por material."""
    for mat in ['fe', 'cr']:
        if datos_totales[mat]:
            heights = [d['height'] for d in datos_totales[mat]]
            altura_promedio = np.mean(heights, axis=0)
            n = altura_promedio.shape[0]

            plt.figure(figsize=(10, 8))
            plt.imshow(altura_promedio, cmap='viridis', interpolation='nearest', 
                      vmin=0, vmax=np.max(altura_promedio))
            plt.colorbar(label='Altura promedio (capas)')
            plt.title(f'Mapa 2D de Altura Promedio - {mat.upper()}')
            plt.xlabel('Coordenada x (sitios)')
            plt.ylabel('Coordenada y (sitios)')
            plt.savefig(os.path.join('./graficos', f'mapa_2d_promedio_{mat}.png'))
            plt.close()

            x = np.arange(n)
            y = np.arange(n)
            X, Y = np.meshgrid(x, y)
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X, Y, altura_promedio, cmap='viridis', linewidth=0, 
                                  antialiased=False, vmin=0, vmax=np.max(altura_promedio))
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
    else:
        print("No se encontraron datos para graficar.")

if __name__ == "__main__":
    main()
