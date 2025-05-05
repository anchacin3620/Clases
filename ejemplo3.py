# Programa para clasificar un delito según su tipo
delito = input("¿Qué tipo de delito se reporta? ").lower()  # Convertir a minúsculas para evitar errores

if delito == "robo":
    print("Enviar patrulla inmediata.")
elif delito == "violencia":
    print("Activar protocolo con defensoría.")
else:
    print("Registrar y clasificar para análisis.")
