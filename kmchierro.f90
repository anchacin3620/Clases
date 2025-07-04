! Programa de Simulación Monte Carlo elaborada por Ángel Chacín Ávila, bajo GNU Fortran y con licencia GPL
PROGRAM simulacion_kmc_pulverizacion
      IMPLICIT NONE
      ! Definición de parámetros constantes
      INTEGER, PARAMETER :: TAM_REJILLA = 300, MAX_PARTICULAS = 100000
      REAL, PARAMETER :: KB = 8.617E-5, ENERGIA_ACTIVACION = 0.7, TEMP_SUBSTRATO = 500.0, TASA_EVENTOS = 100000.0
      ! Declaración de variables
      REAL, DIMENSION(TAM_REJILLA,TAM_REJILLA) :: altura
      INTEGER :: i, j, num_eventos, semilla(8)
      REAL :: tiempo

      ! Validar parámetros físicos para evitar valores no físicos
      IF (ENERGIA_ACTIVACION <= 0.0 .OR. TEMP_SUBSTRATO <= 0.0 .OR. TASA_EVENTOS <= 0.0) THEN
          PRINT *, 'Error: ENERGIA_ACTIVACION, TEMP_SUBSTRATO y TASA_EVENTOS deben ser positivos'
          STOP
      END IF

      ! Inicializar semilla para números aleatorios usando el reloj del sistema
      CALL SYSTEM_CLOCK(i)
      semilla = i + 37 * (/ (j - 1, j = 1, 8) /)
      CALL RANDOM_SEED(PUT=semilla)

      ! Inicializar la rejilla con alturas iniciales aleatorias
      CALL inicializar_rejilla(altura, TAM_REJILLA)

      ! Ejecutar simulación con distribución de Poisson hasta alcanzar MAX_PARTICULAS
      num_eventos = 0
      tiempo = 0.0
      DO WHILE (num_eventos < MAX_PARTICULAS)
          CALL simular_evento(altura, num_eventos, MAX_PARTICULAS, TASA_EVENTOS, tiempo, &
                              ENERGIA_ACTIVACION, TEMP_SUBSTRATO, KB, TAM_REJILLA)
          IF (num_eventos >= MAX_PARTICULAS) EXIT
      END DO

      ! Guardar los resultados en un archivo para visualización con Gnuplot
      CALL guardar_datos(altura, TAM_REJILLA, 'fe_2d_altura.dat', ENERGIA_ACTIVACION, TEMP_SUBSTRATO, TASA_EVENTOS)

      ! Mostrar mensaje de finalización
      PRINT *, 'Simulación 2D para Fe completada. Ver fe_2d_altura.dat'
      STOP

      CONTAINS

      ! Subrutina para inicializar la rejilla con alturas aleatorias
      SUBROUTINE inicializar_rejilla(altura, tam)
          REAL, DIMENSION(tam,tam), INTENT(OUT) :: altura
          INTEGER, INTENT(IN) :: tam
          INTEGER :: i, j
          REAL :: r
          ! Establecer todas las alturas en 0.0
          altura = 0.0
          ! Asignar altura 1.0 a un 5% de los sitios aleatoriamente
          DO i = 1, tam
              DO j = 1, tam
                  CALL RANDOM_NUMBER(r)
                  IF (r < 0.05) altura(i,j) = 1.0
              END DO
          END DO
      END SUBROUTINE inicializar_rejilla

      ! Subrutina para simular un evento (deposición, difusión o pulverización)
      SUBROUTINE simular_evento(altura, num_eventos, max_part, tasa, tiempo, energia_act, temp, kb, tam)
          REAL, DIMENSION(:,:), INTENT(INOUT) :: altura
          INTEGER, INTENT(INOUT) :: num_eventos
          INTEGER, INTENT(IN) :: max_part, tam
          REAL, INTENT(IN) :: tasa, energia_act, temp, kb
          REAL, INTENT(INOUT) :: tiempo
          INTEGER :: x, y, nuevo_x, nuevo_y, direccion
          REAL :: angulo, energia, prob_ads, prob_pulv, r, r_pulv, prob_dif, r_tiempo
          REAL, PARAMETER :: pi = 3.14159265359, MIN_R = 1.0E-10

          ! Generar nuevo evento con distribución de Poisson
          CALL RANDOM_NUMBER(r)
          num_eventos = num_eventos + MAX(1, FLOOR(-LOG(MAX(MIN_R, r)) * tasa))
          IF (num_eventos > max_part) num_eventos = max_part

          ! Generar tiempo entre eventos (distribución exponencial)
          CALL RANDOM_NUMBER(r_tiempo)
          r_tiempo = MAX(MIN_R, r_tiempo)
          tiempo = tiempo + (-LOG(r_tiempo) / tasa)

          ! Simular evento si no se excede el límite de partículas
          IF (num_eventos <= max_part) THEN
              ! Seleccionar sitio aleatorio en la rejilla
              CALL RANDOM_NUMBER(r)
              x = FLOOR(r * (tam-1)) + 1
              CALL RANDOM_NUMBER(r)
              y = FLOOR(r * (tam-1)) + 1
              ! Generar ángulo y energía de la partícula incidente
              CALL RANDOM_NUMBER(r)
              angulo = r * pi / 2.0
              CALL RANDOM_NUMBER(r)
              energia = r * 10.0
              ! Calcular probabilidad de adsorción
              prob_ads = 1.0 - SIN(angulo) * (energia / 10.0)
              CALL RANDOM_NUMBER(r)
              IF (prob_ads > r) THEN
                  ! Deposición: incrementar altura en el sitio
                  altura(x,y) = altura(x,y) + 1.0
                  ! Calcular probabilidad de difusión (Arrhenius)
                  prob_dif = MIN(1.0, EXP(-energia_act / (kb * temp)))
                  CALL RANDOM_NUMBER(r)
                  IF (prob_dif > r) THEN
                      ! Difusión: mover partícula a un sitio vecino
                      CALL RANDOM_NUMBER(r)
                      direccion = FLOOR(r * 4.0)
                      nuevo_x = x
                      nuevo_y = y
                      SELECT CASE (direccion)
                          CASE (0) ! Arriba
                              nuevo_y = y - 1
                              IF (nuevo_y < 1) nuevo_y = tam
                          CASE (1) ! Abajo
                              nuevo_y = y + 1
                              IF (nuevo_y > tam) nuevo_y = 1
                          CASE (2) ! Izquierda
                              nuevo_x = x - 1
                              IF (nuevo_x < 1) nuevo_x = tam
                          CASE (3) ! Derecha
                              nuevo_x = x + 1
                              IF (nuevo_x > tam) nuevo_x = 1
                      END SELECT
                      ! Actualizar alturas para difusión
                      altura(nuevo_x, nuevo_y) = altura(nuevo_x, nuevo_y) + 1.0
                      altura(x, y) = altura(x, y) - 1.0
                  END IF
              ELSE
                  ! Calcular probabilidad de pulverización
                  prob_pulv = SIN(angulo) * (energia / 10.0)
                  CALL RANDOM_NUMBER(r_pulv)
                  IF (prob_pulv > r_pulv .AND. altura(x,y) > 0) THEN
                      ! Pulverización: reducir altura en el sitio
                      altura(x,y) = altura(x,y) - 1.0
                  END IF
              END IF
          END IF
      END SUBROUTINE simular_evento

      ! Subrutina para guardar los datos de altura en un archivo
      SUBROUTINE guardar_datos(altura, tam, nombre_archivo, energia_act, temp, tasa)
          REAL, DIMENSION(tam,tam), INTENT(IN) :: altura
          INTEGER, INTENT(IN) :: tam
          CHARACTER(LEN=*), INTENT(IN) :: nombre_archivo
          REAL, INTENT(IN) :: energia_act, temp, tasa
          INTEGER :: i, j
          ! Abrir archivo para escritura, sobrescribiendo si existe
          OPEN(UNIT=10, FILE=nombre_archivo, STATUS='REPLACE')
          ! Escribir encabezado con parámetros de la simulación
          WRITE(10, '(A)') '# Simulación KMC 2D para pulverización'
          WRITE(10, '(A, I0)') '# Tamaño de la rejilla: ', tam
          WRITE(10, '(A, I0)') '# Número de eventos: ', MAX_PARTICULAS
          WRITE(10, '(A, F8.3)') '# Energía de activación (eV): ', energia_act
          WRITE(10, '(A, F8.1)') '# Temperatura del sustrato (K): ', temp
          WRITE(10, '(A, F8.1)') '# Tasa de eventos: ', tasa
          WRITE(10, '(A)') '# Formato: i j altura(i,j)'
          ! Escribir datos de altura para cada punto de la rejilla
          DO i = 1, tam
              DO j = 1, tam
                  WRITE(10, '(I4, I4, F8.2)') i, j, altura(i,j)
              END DO
          END DO
          ! Cerrar archivo
          CLOSE(10)
      END SUBROUTINE guardar_datos
END PROGRAM simulacion_kmc_pulverizacion
