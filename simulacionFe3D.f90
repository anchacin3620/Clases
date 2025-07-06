! Programa de Simulación Monte Carlo elaborado por Ángel Chacín Ávila, codificado en GNU Fortran y bajo licencia GPL
PROGRAM simulacion_montecarlo_pulverizacion_3d_fe
      IMPLICIT NONE
      INTEGER, PARAMETER :: TAM_REJILLA = 300, ALTURA_MAX = 100, MAX_PARTICULAS = 3000000
      REAL, PARAMETER :: KB = 8.617E-5, ENERGIA_DIFUSION = 0.7, TEMP_SUBSTRATO = 500.0, TEMP_PLASMA = 1000.0
      REAL, PARAMETER :: ADHESION = 0.85, TASA_EVENTOS = 3000000.0
      INTEGER, DIMENSION(TAM_REJILLA,TAM_REJILLA,ALTURA_MAX) :: material
      INTEGER :: i, j, k, num_eventos, semilla(8)
      REAL :: tiempo

      ! Validar parámetros físicos
      IF (ENERGIA_DIFUSION <= 0.0 .OR. TEMP_SUBSTRATO <= 0.0 .OR. TEMP_PLASMA <= 0.0 .OR. TASA_EVENTOS <= 0.0 .OR. ADHESION <= 0.0) THEN
          PRINT *, 'Error: ENERGIA_DIFUSION, TEMP_SUBSTRATO, TEMP_PLASMA, TASA_EVENTOS y ADHESION deben ser positivos'
          STOP
      END IF

      ! Inicializar semilla para números aleatorios (portable)
      CALL SYSTEM_CLOCK(i)
      semilla = i + 37 * (/ (j - 1, j = 1, 8) /)
      CALL RANDOM_SEED(PUT=semilla)

      ! Inicializar rejilla
      CALL inicializar_rejilla(material, TAM_REJILLA, ALTURA_MAX)

      ! Simulación con distribución de Poisson
      num_eventos = 0
      tiempo = 0.0
      DO WHILE (num_eventos < MAX_PARTICULAS)
          CALL simular_evento(material, num_eventos, MAX_PARTICULAS, TASA_EVENTOS, tiempo, ENERGIA_DIFUSION, TEMP_SUBSTRATO, KB, ADHESION, TEMP_PLASMA, TAM_REJILLA, ALTURA_MAX)
          IF (num_eventos >= MAX_PARTICULAS) EXIT
      END DO

      ! Guardar resultados para Gnuplot
      CALL guardar_datos(material, TAM_REJILLA, ALTURA_MAX, 'fe_3d_altura.dat', ENERGIA_DIFUSION, TEMP_SUBSTRATO, TASA_EVENTOS, ADHESION)

      PRINT *, 'Simulación 3D para Fe completada. Ver fe_3d_altura.dat'
      STOP

      CONTAINS

      SUBROUTINE inicializar_rejilla(material, tam_rejilla, altura_max)
          INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(OUT) :: material
          INTEGER, INTENT(IN) :: tam_rejilla, altura_max
          INTEGER :: i, j
          material = 0
          DO i = 1, tam_rejilla
              DO j = 1, tam_rejilla
                  material(i,j,1) = 3
                  material(i,j,2) = 3
              END DO
          END DO
      END SUBROUTINE inicializar_rejilla

      SUBROUTINE simular_evento(material, num_eventos, max_part, tasa_eventos, tiempo, energia_difusion, temp_substrato, kb, adhesion, temp_plasma, tam_rejilla, altura_max)
          INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: material
          INTEGER, INTENT(INOUT) :: num_eventos
          INTEGER, INTENT(IN) :: max_part, tam_rejilla, altura_max
          REAL, INTENT(IN) :: tasa_eventos, energia_difusion, temp_substrato, kb, adhesion, temp_plasma
          REAL, INTENT(INOUT) :: tiempo
          INTEGER :: x, y, z_max, nuevo_x, nuevo_y, direccion
          REAL :: energia, prob_adhesion, r, prob_difusion, r_tiempo
          REAL, PARAMETER :: pi = 3.14159265359, MIN_R = 1.0E-10

          CALL RANDOM_NUMBER(r)
          num_eventos = num_eventos + MAX(1, FLOOR(-LOG(MAX(MIN_R, r)) * tasa_eventos))
          IF (num_eventos > max_part) num_eventos = max_part

          CALL RANDOM_NUMBER(r_tiempo)
          r_tiempo = MAX(MIN_R, r_tiempo)
          tiempo = tiempo + (-LOG(r_tiempo) / tasa_eventos)

          IF (num_eventos <= max_part) THEN
              CALL RANDOM_NUMBER(r)
              x = FLOOR(r * (tam_rejilla-1)) + 1
              CALL RANDOM_NUMBER(r)
              y = FLOOR(r * (tam_rejilla-1)) + 1
              z_max = 2
              DO i = MAX(1,x-1), MIN(tam_rejilla,x+1)
                  DO j = MAX(1,y-1), MIN(tam_rejilla,y+1)
                      z_max = MAX(z_max, MAXVAL(material(i,j,1:altura_max)))
                  END DO
              END DO
              CALL RANDOM_NUMBER(r)
              energia = -KB * temp_plasma * LOG(MAX(MIN_R, r))
              prob_adhesion = adhesion * EXP(-energia / (KB * temp_plasma))
              CALL RANDOM_NUMBER(r)
              IF (prob_adhesion > r) THEN
                  material(x,y,z_max+1) = 2
                  prob_difusion = MIN(1.0, EXP(-energia_difusion / (kb * temp_substrato)))
                  CALL RANDOM_NUMBER(r)
                  IF (prob_difusion > r .AND. z_max > 2) THEN
                      CALL RANDOM_NUMBER(r)
                      direccion = FLOOR(r * 4.0)
                      nuevo_x = x
                      nuevo_y = y
                      SELECT CASE (direccion)
                          CASE (0) ! Arriba
                              nuevo_y = y - 1
                              IF (nuevo_y < 1) nuevo_y = tam_rejilla
                          CASE (1) ! Abajo
                              nuevo_y = y + 1
                              IF (nuevo_y > tam_rejilla) nuevo_y = 1
                          CASE (2) ! Izquierda
                              nuevo_x = x - 1
                              IF (nuevo_x < 1) nuevo_x = tam_rejilla
                          CASE (3) ! Derecha
                              nuevo_x = x + 1
                              IF (nuevo_x > tam_rejilla) nuevo_x = 1
                      END SELECT
                      material(nuevo_x, nuevo_y, z_max+1) = 2
                      material(x, y, z_max+1) = 0
                  END IF
              END IF
          END IF
      END SUBROUTINE simular_evento

      SUBROUTINE guardar_datos(material, tam_rejilla, altura_max, nombre_archivo, energia_difusion, temp_substrato, tasa_eventos, adhesion)
          INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(IN) :: material
          INTEGER, INTENT(IN) :: tam_rejilla, altura_max
          CHARACTER(LEN=*), INTENT(IN) :: nombre_archivo
          REAL, INTENT(IN) :: energia_difusion, temp_substrato, tasa_eventos, adhesion
          INTEGER :: i, j, z_max
          OPEN(UNIT=10, FILE=nombre_archivo, STATUS='REPLACE')
          WRITE(10, '(A)') '# Simulación KMC 3D para pulverización (Fe)'
          WRITE(10, '(A, I0)') '# Tamaño de la rejilla (TAM_REJILLA): ', tam_rejilla
          WRITE(10, '(A, I0)') '# Altura máxima (ALTURA_MAX): ', altura_max
          WRITE(10, '(A, I0)') '# Número de eventos (MAX_PARTICULAS): ', MAX_PARTICULAS
          WRITE(10, '(A, F8.3)') '# Energía de difusión (eV): ', energia_difusion
          WRITE(10, '(A, F8.1)') '# Temperatura del sustrato (K): ', temp_substrato
          WRITE(10, '(A, F8.1)') '# Tasa de eventos: ', tasa_eventos
          WRITE(10, '(A, F8.3)') '# Probabilidad de adhesión (ADHESION): ', adhesion
          WRITE(10, '(A)') '# Formato: i j altura(i,j)'
          DO i = 1, tam_rejilla
              DO j = 1, tam_rejilla
                  z_max = MAXVAL(material(i,j,1:altura_max))
                  WRITE(10, '(I4, I4, I4)') i, j, z_max
              END DO
          END DO
          CLOSE(10)
      END SUBROUTINE guardar_datos
END PROGRAM simulacion_montecarlo_pulverizacion_3d_fe
