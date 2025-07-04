import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Cargar datos del archivo (ajusta la ruta según tu ubicación)
data = np.loadtxt('fe_2d_altura.dat')  # Asume formato: i j altura(i,j) en tres columnas
n = int(np.sqrt(len(data)))  # Tamaño de la rejilla (300x300)
height = data[:, 2].reshape(n, n)  # Reorganiza las alturas en una matriz 300x300

# Calcular estadísticas
rugosidad = np.std(height)  # Desviación estándar (rugosidad)
espesor_promedio = np.mean(height)  # Media de las alturas
max_altura = np.max(height)  # Altura máxima
min_altura = np.min(height)  # Altura mínima

print(f"Rugosidad (desviación estándar): {rugosidad:.4f}")
print(f"Espesor promedio: {espesor_promedio:.4f}")
print(f"Altura máxima: {max_altura:.4f}")
print(f"Altura mínima: {min_altura:.4f}")

# Visualización 2D
plt.figure(figsize=(10, 8))
plt.imshow(height, cmap='viridis', interpolation='nearest', vmin=min_altura, vmax=max_altura)
plt.colorbar(label='Altura')
plt.title('Mapa 2D de Altura - Depósito de 100000 Partículas Fe')
plt.xlabel('Coordenada x (sitios)')
plt.ylabel('Coordenada y (sitios)')
plt.show()

# Visualización 3D
x = np.arange(n)
y = np.arange(n)
X, Y = np.meshgrid(x, y)
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, height, cmap='viridis', linewidth=0, antialiased=False)
ax.set_zlim(min_altura, max_altura)
ax.set_xlabel('Coordenada x (sitios)')
ax.set_ylabel('Coordenada y (sitios)')
ax.set_zlabel('Altura')
ax.view_init(30, 45)  # Ángulo de vista ajustado
plt.title('Mapa 3D de Altura - Depósito de 100000 Partículas Fe')
plt.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

# Identificar islas (regiones con altura > umbral, e.g., 0.6)
umbral = 0.6
islas = height > umbral
num_islas, etiquetas = np.unique(islas, return_counts=True)
cobertura = np.sum(islas) / (n * n) * 100  # Porcentaje de cobertura
print(f"Cobertura de islas (> {umbral}): {cobertura:.2f}%")
