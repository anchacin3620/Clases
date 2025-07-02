      PROGRAM MonteCarloSputtering3D_Cr
      IMPLICIT NONE
      INTEGER, PARAMETER :: L = 300, H = 100, MAX_PART = 20000
      REAL, PARAMETER :: KB = 8.617E-5, ED = 0.6, T_SUBS = 500.0, T_PLASMA = 1000.0
      REAL, PARAMETER :: STICKING = 0.9, LAMBDA = 20000.0
      INTEGER, DIMENSION(L,L,H) :: material
      INTEGER :: i, j, k, n_events, seed(8)
      REAL :: time

      ! Validar parámetros físicos
      IF (ED <= 0.0 .OR. T_SUBS <= 0.0 .OR. T_PLASMA <= 0.0 .OR. LAMBDA <= 0.0 .OR. STICKING <= 0.0) THEN
          PRINT *, 'Error: ED, T_SUBS, T_PLASMA, LAMBDA y STICKING deben ser positivos'
          STOP
      END IF

      ! Inicializar semilla para números aleatorios (portable)
      CALL SYSTEM_CLOCK(i)
      seed = i + 37 * (/ (j - 1, j = 1, 8) /)
      CALL RANDOM_SEED(PUT=seed)

      ! Inicializar rejilla
      CALL initialize_grid(material, L, H)

      ! Simulación con distribución de Poisson
      n_events = 0
      time = 0.0
      DO WHILE (n_events < MAX_PART)
          CALL simulate_event(material, n_events, MAX_PART, LAMBDA, time, ED, T_SUBS, KB, STICKING, T_PLASMA, L, H)
          IF (n_events >= MAX_PART) EXIT
      END DO

      ! Guardar resultados para Gnuplot
      CALL save_data(material, L, H, 'cr_3d_height.dat', ED, T_SUBS, LAMBDA, STICKING)

      PRINT *, 'Simulación 3D para Cr completada. Ver cr_3d_height.dat'
      STOP
      CONTAINS

      SUBROUTINE initialize_grid(material, L, H)
          INTEGER, DIMENSION(L,L,H), INTENT(OUT) :: material
          INTEGER, INTENT(IN) :: L, H
          INTEGER :: i, j
          material = 0
          DO i = 1, L
              DO j = 1, L
                  material(i,j,1) = 3
                  material(i,j,2) = 3
              END DO
          END DO
      END SUBROUTINE initialize_grid

      SUBROUTINE simulate_event(material, n_events, max_part, lambda, time, ed, t_subs, kb, sticking, t_plasma, L, H)
          INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: material
          INTEGER, INTENT(INOUT) :: n_events
          INTEGER, INTENT(IN) :: max_part, L, H
          REAL, INTENT(IN) :: lambda, ed, t_subs, kb, sticking, t_plasma
          REAL, INTENT(INOUT) :: time
          INTEGER :: x, y, z_max, new_x, new_y, direction
          REAL :: E, P_stick, r, P_diff, r_time
          REAL, PARAMETER :: pi = 3.14159265359, MIN_R = 1.0E-10

          CALL RANDOM_NUMBER(r)
          n_events = n_events + MAX(1, FLOOR(-LOG(MAX(MIN_R, r)) * lambda))
          IF (n_events > max_part) n_events = max_part

          CALL RANDOM_NUMBER(r_time)
          r_time = MAX(MIN_R, r_time)
          time = time + (-LOG(r_time) / lambda)

          IF (n_events <= max_part) THEN
              CALL RANDOM_NUMBER(r)
              x = FLOOR(r * (L-1)) + 1
              CALL RANDOM_NUMBER(r)
              y = FLOOR(r * (L-1)) + 1
              z_max = 2
              DO i = MAX(1,x-1), MIN(L,x+1)
                  DO j = MAX(1,y-1), MIN(L,y+1)
                      z_max = MAX(z_max, MAXVAL(material(i,j,1:H)))
                  END DO
              END DO
              CALL RANDOM_NUMBER(r)
              E = -KB * t_plasma * LOG(MAX(MIN_R, r))
              P_stick = sticking * EXP(-E / (KB * t_plasma))
              CALL RANDOM_NUMBER(r)
              IF (P_stick > r) THEN
                  material(x,y,z_max+1) = 1
                  P_diff = MIN(1.0, EXP(-ed / (kb * t_subs)))
                  CALL RANDOM_NUMBER(r)
                  IF (P_diff > r .AND. z_max > 2) THEN
                      CALL RANDOM_NUMBER(r)
                      direction = FLOOR(r * 4.0)
                      new_x = x
                      new_y = y
                      SELECT CASE (direction)
                          CASE (0) ! Arriba
                              new_y = y - 1
                              IF (new_y < 1) new_y = L
                          CASE (1) ! Abajo
                              new_y = y + 1
                              IF (new_y > L) new_y = 1
                          CASE (2) ! Izquierda
                              new_x = x - 1
                              IF (new_x < 1) new_x = L
                          CASE (3) ! Derecha
                              new_x = x + 1
                              IF (new_x > L) new_x = 1
                      END SELECT
                      material(new_x, new_y, z_max+1) = 1
                      material(x, y, z_max+1) = 0
                  END IF
              END IF
          END IF
      END SUBROUTINE simulate_event

      SUBROUTINE save_data(material, L, H, filename, ed, t_subs, lambda, sticking)
          INTEGER, DIMENSION(L,L,H), INTENT(IN) :: material
          INTEGER, INTENT(IN) :: L, H
          CHARACTER(LEN=*), INTENT(IN) :: filename
          REAL, INTENT(IN) :: ed, t_subs, lambda, sticking
          INTEGER :: i, j, z_max
          OPEN(UNIT=10, FILE=filename, STATUS='REPLACE')
          WRITE(10, '(A)') '# Simulación KMC 3D para sputtering (Cr)'
          WRITE(10, '(A, I0)') '# Tamaño de la rejilla (L): ', L
          WRITE(10, '(A, I0)') '# Altura (H): ', H
          WRITE(10, '(A, I0)') '# Número de eventos (MAX_PART): ', MAX_PART
          WRITE(10, '(A, F8.3)') '# Energía de activación (ED, eV): ', ed
          WRITE(10, '(A, F8.1)') '# Temperatura del sustrato (T_SUBS, K): ', t_subs
          WRITE(10, '(A, F8.1)') '# Tasa de eventos (LAMBDA): ', lambda
          WRITE(10, '(A, F8.3)') '# Probabilidad de sticking (STICKING): ', sticking
          WRITE(10, '(A)') '# Formato: i j height(i,j)'
          DO i = 1, L
              DO j = 1, L
                  z_max = MAXVAL(material(i,j,1:H))
                  WRITE(10, '(I4, I4, I4)') i, j, z_max
              END DO
          END DO
          CLOSE(10)
      END SUBROUTINE save_data
      END PROGRAM MonteCarloSputtering3D_Cr
