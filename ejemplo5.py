# Sistema Integrado de Seguridad Ciudadana
print("=== Sistema de Seguridad Ciudadana ===")
incidentes = []  # Lista para almacenar los incidentes

# Revisión de 3 cámaras
for i in range(3):
    print(f"\nRevisando cámara de seguridad número {i+1}")
    delito = input("¿Qué tipo de delito se reporta? ").lower()
    incidentes.append(delito)  # Guardar en lista
    
    # Clasificación
    if delito == "robo":
        print("-> Enviar patrulla inmediata.")
    elif delito == "violencia":
        print("-> Activar protocolo con defensoría.")
    else:
        print("-> Registrar y clasificar para análisis.")

# Reporte de incidentes almacenados
print("\n=== Resumen de incidentes ===")
print("Total de incidentes registrados:", len(incidentes))
print("Primer incidente reportado:", incidentes[0])
print("Todos los incidentes:", incidentes)
