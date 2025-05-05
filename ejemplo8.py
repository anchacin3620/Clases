# Sistema Integrado de Seguridad Ciudadana - Versión Ampliada
print("=== SISTEMA DE MONITOREO URBANO INTELIGENTE ===")
incidentes = []
zonas = []
alertas = []
codigos = []

# Catálogo de respuestas
protocolos = {
    ("Centro", "secuestro"): ["NEGRA", "Movilizar BAE y fiscalía", "P-001"],
    ("Transporte", "robo"): ["AMARILLA", "Bloquear salidas y revisar CCTV", "P-002"],
    ("Residencial", "violencia"): ["MORADA", "Protocolo de violencia doméstica", "P-003"],
    ("Comercial", "vandalismo"): ["VERDE", "Patrulla rápida + reporte a seguros", "P-004"],
    ("Escolar", "amenaza"): ["ROJA", "Evacuación preventiva + EOD", "P-005"],  # Nueva combinación 1
    ("Hospital", "asalto"): ["AZUL", "Bloqueo médico + intervención táctica", "P-006"],  # Nueva combinación 2
    ("Turística", "fraude"): ["DORADA", "Intervención financiera + interpol", "P-007"],  # Nueva combinación 3
    ("default"): ["GRIS", "Registrar para análisis", "P-000"]
}

for i in range(3):
    print(f"\n?? Cámara {i+1}/3")
    zona = input("?? Zona (Centro/Transporte/Residencial/Comercial/Escolar/Hospital/Turística): ").capitalize()
    delito = input("?? Incidente (robo/violencia/secuestro/vandalismo/amenaza/asalto/fraude): ").lower()
    
    # Búsqueda de protocolo
    accion = protocolos.get((zona, delito), protocolos["default"])
    
    # Registro
    incidentes.append(delito)
    zonas.append(zona)
    alertas.append(accion[0])
    codigos.append(accion[2])
    
    print(f"\n?? ALERTA {accion[0]}: {accion[1]} (Código: {accion[2]})")

# Reporte mejorado
print("\n" + "="*60)
print(f"{'?? INFORME FINAL':^60}")
print("="*60)
print(f"{'Cámara':<8}{'Zona':<15}{'Incidente':<15}{'Alerta':<10}{'Protocolo'}")
print("-"*60)
for i in range(3):
    print(f"{i+1:<8}{zonas[i]:<15}{incidentes[i]:<15}{alertas[i]:<10}{codigos[i]}")
print("="*60)
