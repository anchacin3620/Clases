      PROGRAM kmc_sputtering_simulation
      IMPLICIT NONE
      INTEGER, PARAMETER :: L = 300, MAX_PART = 1000
      REAL, PARAMETER :: KB = 8.617E-5, ED = 0.6, T_SUBS = 500.0, LAMBDA = 1000.0
      REAL, DIMENSION(L,L) :: height
      INTEGER :: i, j, n_events, seed(8)
      REAL :: time

      ! Validar parámetros físicos
      IF (ED <= 0.0 .OR. T_SUBS <= 0.0 .OR. LAMBDA <= 0.0) THEN
          PRINT *, 'Error: ED, T_SUBS y LAMBDA deben ser positivos'
          STOP
      END IF

      ! Inicializar semilla para números aleatorios (portable)
      CALL SYSTEM_CLOCK(i)
      seed = i + 37 * (/ (j - 1, j = 1, 8) /)
      CALL RANDOM_SEED(PUT=seed)

      ! Inicializar rejilla
      CALL initialize_grid(height, L)

      ! Simulación con distribución de Poisson
      n_events = 0
      time = 0.0
      DO WHILE (n_events < MAX_PART)
          CALL simulate_event(height, n_events, MAX_PART, LAMBDA, time, ED, T_SUBS, KB, L)
          IF (n_events >= MAX_PART) EXIT
      END DO

      ! Guardar resultados para Gnuplot
      CALL save_data(height, L, 'cr_2d_height.dat', ED, T_SUBS, LAMBDA)

      PRINT *, 'Simulación 2D para Cr completada. Ver cr_2d_height.dat'
      STOP
      CONTAINS

      SUBROUTINE initialize_grid(height, L)
          REAL, DIMENSION(L,L), INTENT(OUT) :: height
          INTEGER, INTENT(IN) :: L
          INTEGER :: i, j
          REAL :: r
          height = 0.0
          DO i = 1, L
              DO j = 1, L
                  CALL RANDOM_NUMBER(r)
                  IF (r < 0.05) height(i,j) = 1.0
              END DO
          END DO
      END SUBROUTINE initialize_grid

      SUBROUTINE simulate_event(height, n_events, max_part, lambda, time, ed, t_subs, kb, L)
          REAL, DIMENSION(:,:), INTENT(INOUT) :: height
          INTEGER, INTENT(INOUT) :: n_events
          INTEGER, INTENT(IN) :: max_part, L
          REAL, INTENT(IN) :: lambda, ed, t_subs, kb
          REAL, INTENT(INOUT) :: time
          INTEGER :: x, y, new_x, new_y, direction
          REAL :: theta, E, P_ads, P_sput, r, r_sput, P_diff, r_time
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
              CALL RANDOM_NUMBER(r)
              theta = r * pi / 2.0
              CALL RANDOM_NUMBER(r)
              E = r * 10.0
              P_ads = 1.0 - SIN(theta) * (E / 10.0)
              CALL RANDOM_NUMBER(r)
              IF (P_ads > r) THEN
                  height(x,y) = height(x,y) + 1.0
                  P_diff = MIN(1.0, EXP(-ed / (kb * t_subs)))
                  CALL RANDOM_NUMBER(r)
                  IF (P_diff > r) THEN
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
                      height(new_x, new_y) = height(new_x, new_y) + 1.0
                      height(x, y) = height(x, y) - 1.0
                  END IF
              ELSE
                  P_sput = SIN(theta) * (E / 10.0)
                  CALL RANDOM_NUMBER(r_sput)
                  IF (P_sput > r_sput .AND. height(x,y) > 0) THEN
                      height(x,y) = height(x,y) - 1.0
                  END IF
              END IF
          END IF
      END SUBROUTINE simulate_event

      SUBROUTINE save_data(height, L, filename, ed, t_subs, lambda)
          REAL, DIMENSION(L,L), INTENT(IN) :: height
          INTEGER, INTENT(IN) :: L
          CHARACTER(LEN=*), INTENT(IN) :: filename
          REAL, INTENT(IN) :: ed, t_subs, lambda
          INTEGER :: i, j
          OPEN(UNIT=10, FILE=filename, STATUS='REPLACE')
          WRITE(10, '(A)') '# Simulación KMC 2D para sputtering'
          WRITE(10, '(A, I0)') '# Tamaño de la rejilla (L): ', L
          WRITE(10, '(A, I0)') '# Número de eventos (MAX_PART): ', MAX_PART
          WRITE(10, '(A, F8.3)') '# Energía de activación (ED, eV): ', ed
          WRITE(10, '(A, F8.1)') '# Temperatura del sustrato (T_SUBS, K): ', t_subs
          WRITE(10, '(A, F8.1)') '# Tasa de eventos (LAMBDA): ', lambda
          WRITE(10, '(A)') '# Formato: i j height(i,j)'
          DO i = 1, L
              DO j = 1, L
                  WRITE(10, '(I4, I4, F8.2)') i, j, height(i,j)
              END DO
          END DO
          CLOSE(10)
      END SUBROUTINE save_data
      END PROGRAM kmc_sputtering_simulation
