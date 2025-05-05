# Sistema Integrado de Seguridad Ciudadana con Priorización
print("=== Sistema de Seguridad Ciudadana Avanzado ===")
incidentes = []  # Lista para almacenar los incidentes
zonas = []       # Lista para almacenar las zonas

# Revisión de 3 cámaras
for i in range(3):
    print(f"\nRevisando cámara de seguridad número {i+1}")
    zona = input("Zona del incidente (Centro/Periferia/Otros): ").capitalize()
    delito = input("¿Qué tipo de delito se reporta? ").lower()
    
    # Almacenamiento de datos
    incidentes.append(delito)
    zonas.append(zona)
    
    # Sistema de priorización
    if zona == "Centro" and delito == "violencia":
        print("¡Alerta roja! Enviar unidades especiales.")
    elif zona == "Centro" and delito == "robo":
        print("¡Alerta naranja! Enviar patrulla reforzada.")
    elif delito == "violencia":
        print("Activar protocolo con defensoría.")
    elif delito == "robo":
        print("Enviar patrulla estándar.")
    else:
        print("Registrar y clasificar para análisis.")

# Reporte de incidentes almacenados
print("\n=== Resumen de incidentes ===")
print(f"Total de incidentes registrados: {len(incidentes)}")
print("Detalle por zona y tipo:")
for i in range(len(incidentes)):
    print(f"Cámara {i+1}: Zona {zonas[i]} - {incidentes[i]}")
