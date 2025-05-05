# Programa que simula revisión de cámaras y clasificación de delitos
print("=== Sistema de Seguridad Ciudadana ===")

# Simular revisión de 3 cámaras
for i in range(3):
    print("\nRevisando cámara de seguridad número", i + 1)
    delito = input("¿Qué tipo de delito se reporta? ").lower()
    
    if delito == "robo":
        print("-> Enviar patrulla inmediata.")
    elif delito == "violencia":
        print("-> Activar protocolo con defensoría.")
    else:
        print("-> Registrar y clasificar para análisis.")

print("\nProceso de revisión completado.")
