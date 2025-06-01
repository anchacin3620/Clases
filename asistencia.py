import datetime
import csv
import os
from pathlib import Path

class SistemaAsistencia:
    def __init__(self):
        self.registros = []
        self.ruta_archivo = "asistencia_unes_zulia.csv"

    def validar_cedula(self, cedula):
        """Valida el formato de la cédula venezolana"""
        try:
            return cedula.isdigit() and 6 <= len(cedula) <= 8
        except:
            return False

    def registrar_asistencia(self):
        """Registra la asistencia de un trabajador"""
        print("\n=== Registro de Asistencia ===")
        nombre = input("Ingrese nombre: ").strip().title()
        apellido = input("Ingrese apellido: ").strip().title()
        
        while True:
            cedula = input("Ingrese número de cédula (sin puntos ni letras): ").strip()
            if self.validar_cedula(cedula):
                break
            print("Cédula inválida. Debe contener solo números (6-8 dígitos).")

        tipo_personal = input("Tipo de personal (Obrero/Administrativo/Docente): ").strip().title()
        while tipo_personal not in ["Obrero", "Administrativo", "Docente"]:
            print("Tipo de personal inválido.")
            tipo_personal = input("Tipo de personal (Obrero/Administrativo/Docente): ").strip().title()

        fecha = datetime.datetime.now().strftime("%Y-%m-%d")
        hora = datetime.datetime.now().strftime("%H:%M:%S")
        estado = "Presente"

        registro = {
            "Fecha": fecha,
            "Hora": hora,
            "Nombre": nombre,
            "Apellido": apellido,
            "Cedula": cedula,
            "Tipo_Personal": tipo_personal,
            "Estado": estado
        }
        self.registros.append(registro)
        self.guardar_registro(registro)
        print("\nAsistencia registrada exitosamente!")

    def guardar_registro(self, registro):
        """Guarda el registro en un archivo CSV"""
        archivo_existe = Path(self.ruta_archivo).is_file()
        
        with open(self.ruta_archivo, mode='a', newline='', encoding='utf-8') as file:
            writer = csv.DictWriter(file, 
                                  fieldnames=["Fecha", "Hora", "Nombre", "Apellido", 
                                            "Cedula", "Tipo_Personal", "Estado"])
            if not archivo_existe:
                writer.writeheader()
            writer.writerow(registro)

    def generar_reporte(self):
        """Genera un reporte de asistencia por cédula"""
        cedula = input("\nIngrese la cédula del trabajador: ").strip()
        if not self.validar_cedula(cedula):
            print("Cédula inválida.")
            return

        print(f"\n=== Reporte de Asistencia para Cédula: {cedula} ===")
        asistencias = 0
        inasistencias = 0
        registros_encontrados = []

        try:
            with open(self.ruta_archivo, mode='r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row["Cedula"] == cedula:
                        registros_encontrados.append(row)
                        if row["Estado"] == "Presente":
                            asistencias += 1
                        else:
                            inasistencias += 1

            if not registros_encontrados:
                print("No se encontraron registros para esta cédula.")
                return

            # Mostrar información del trabajador
            print(f"Nombre: {registros_encontrados[0]['Nombre']} {registros_encontrados[0]['Apellido']}")
            print(f"Tipo de Personal: {registros_encontrados[0]['Tipo_Personal']}")
            print(f"\nTotal Asistencias: {asistencias}")
            print(f"Total Inasistencias: {inasistencias}")
            print("\nDetalles de registros:")
            for registro in registros_encontrados:
                print(f"Fecha: {registro['Fecha']} | Hora: {registro['Hora']} | Estado: {registro['Estado']}")

        except FileNotFoundError:
            print("No existe un archivo de registros.")

    def registrar_inasistencia(self):
        """Registra una inasistencia manualmente"""
        print("\n=== Registro de Inasistencia ===")
        nombre = input("Ingrese nombre: ").strip().title()
        apellido = input("Ingrese apellido: ").strip().title()
        
        while True:
            cedula = input("Ingrese número de cédula: ").strip()
            if self.validar_cedula(cedula):
                break
            print("Cédula inválida. Debe contener solo números (6-8 dígitos).")

        tipo_personal = input("Tipo de personal (Obrero/Administrativo/Docente): ").strip().title()
        while tipo_personal not in ["Obrero", "Administrativo", "Docente"]:
            print("Tipo de personal inválido.")
            tipo_personal = input("Tipo de personal (Obrero/Administrativo/Docente): ").strip().title()

        fecha = datetime.datetime.now().strftime("%Y-%m-%d")
        hora = datetime.datetime.now().strftime("%H:%M:%S")
        estado = "Ausente"

        registro = {
            "Fecha": fecha,
            "Hora": hora,
            "Nombre": nombre,
            "Apellido": apellido,
            "Cedula": cedula,
            "Tipo_Personal": tipo_personal,
            "Estado": estado
        }
        self.guardar_registro(registro)
        print("\nInasistencia registrada exitosamente!")

    def menu(self):
        """Muestra el menú principal y gestiona las opciones"""
        while True:
            print("\n=== Sistema de Gestión de Asistencia UNES Zulia ===")
            print("1. Registrar Asistencia")
            print("2. Registrar Inasistencia")
            print("3. Generar Reporte de Asistencia")
            print("4. Salir")
            
            opcion = input("\nSeleccione una opción (1-4): ").strip()
            
            if opcion == "1":
                self.registrar_asistencia()
            elif opcion == "2":
                self.registrar_inasistencia()
            elif opcion == "3":
                self.generar_reporte()
            elif opcion == "4":
                print("¡Gracias por usar el sistema!")
                break
            else:
                print("Opción inválida. Por favor, seleccione una opción válida.")

if __name__ == "__main__":
    sistema = SistemaAsistencia()
    sistema.menu()
