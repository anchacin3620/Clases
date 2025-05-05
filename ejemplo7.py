# Sistema Integrado de Seguridad Ciudadana con Priorización Ampliada
print("=== Sistema de Seguridad Ciudadana Avanzado ===")
incidentes = []  # Lista para tipos de delito
zonas = []       # Lista para zonas
alertas = []     # Lista para nivel de alerta

# Revisión de 3 cámaras
for i in range(3):
    print(f"\nRevisando cámara de seguridad número {i+1}")
    zona = input("Zona (Centro/Periferia/Residencial/Comercial/Transporte): ").capitalize()
    delito = input("Tipo de delito (robo/violencia/vandalismo/fraude/secuestro): ").lower()
    
    # Almacenamiento de datos
    incidentes.append(delito)
    zonas.append(zona)
    
    # Sistema de priorización ampliado
    if zona == "Centro" and delito == "secuestro":
        accion = "¡ALERTA NEGRA! Movilizar BAE y fiscalía."
        alertas.append("NEGRA")
    elif zona == "Transporte" and delito == "robo":
        accion = "¡ALERTA AMARILLA! Bloquear salidas y revisar CCTV."
        alertas.append("AMARILLA")
    elif zona == "Residencial" and delito == "violencia":
        accion = "¡ALERTA MORADA! Protocolo de violencia doméstica."
        alertas.append("MORADA")
    elif zona == "Comercial" and delito == "vandalismo":
        accion = "¡ALERTA VERDE! Patrulla rápida y reporte a seguros."
        alertas.append("VERDE")
    elif delito == "fraude":
        accion = "Derivar a unidad de ciberdelitos."
        alertas.append("BLANCA")
    else:
        accion = "Registrar para análisis posterior."
        alertas.append("GRIS")
    
    print(f"-> {accion}")

# Reporte consolidado
print("\n=== Resumen de Incidentes ===")
print(f"{'Cámara':<7} | {'Zona':<12} | {'Delito':<12} | {'Alerta':<10} | Acción")
print("-"*60)
for i in range(len(incidentes)):
    print(f"{i+1:<7} | {zonas[i]:<12} | {incidentes[i]:<12} | {alertas[i]:<10} | ", end="")
    if alertas[i] == "NEGRA":
        print("SWAT + Fiscalía")
    elif alertas[i] == "AMARILLA":
        print("Bloqueo salidas")
    elif alertas[i] == "MORADA":
        print("Violencia doméstica")
    elif alertas[i] == "VERDE":
        print("Patrulla rápida")
    else:
        print("Análisis estándar")
