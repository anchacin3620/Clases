      PROGRAM kmc_sputtering_simulation
      IMPLICIT NONE
      INTEGER, PARAMETER :: L = 300, MAX_PART = 1000
      REAL, PARAMETER :: KB = 8.617E-5, ED = 0.7, T_SUBS = 500.0, LAMBDA = 1000.0
      REAL, DIMENSION(L,L) :: height
      INTEGER :: i, j, n_events, seed
      REAL :: time

      ! Inicializar semilla para números aleatorios
      CALL SYSTEM_CLOCK(seed)
      CALL RANDOM_SEED(seed)

      ! Inicializar rejilla
      CALL initialize_grid(height, L)

      ! Simulación con distribución de Poisson
      n_events = 0
      time = 0.0
      DO WHILE (n_events < MAX_PART)
          CALL simulate_event(height, n_events, MAX_PART, LAMBDA, time, ED, T_SUBS, KB)
          IF (n_events >= MAX_PART) EXIT
      END DO

      ! Guardar resultados para Gnuplot
      CALL save_data(height, L, 'fe_2d_height.dat')

      PRINT *, 'Simulación 2D para Fe completada. Ver fe_2d_height.dat'
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
                  IF (r < 0.05) height(i,j) = 1.0 ! 5% de defectos
              END DO
          END DO
      END SUBROUTINE initialize_grid

      SUBROUTINE simulate_event(height, n_events, max_part, lambda, time, ed, t_subs, kb)
          REAL, DIMENSION(:,:), INTENT(INOUT) :: height
          INTEGER, INTENT(INOUT) :: n_events
          INTEGER, INTENT(IN) :: max_part
          REAL, INTENT(IN) :: lambda, time, ed, t_subs, kb
          INTEGER :: x, y, i, j
          REAL :: theta, E, P_ads, P_sput, r, r_sput, P_diff, r_time
          REAL, PARAMETER :: pi = 3.14159265359

          ! Generar nuevo evento con Poisson
          CALL RANDOM_NUMBER(r)
          n_events = n_events + FLOOR(-LOG(r) * lambda)
          IF (n_events > max_part) n_events = max_part

          ! Generar tiempo entre eventos (exponencial)
          CALL RANDOM_NUMBER(r_time)
          time = time - LOG(r_time) / lambda

          ! Simular evento si no se excede el límite
          IF (n_events <= max_part) THEN
              CALL RANDOM_NUMBER(r)
              x = MIN(MAX(1, INT(r * L) + 1), L)
              CALL RANDOM_NUMBER(r)
              y = MIN(MAX(1, INT(r * L) + 1), L)
              CALL RANDOM_NUMBER(r)
              theta = r * pi / 2.0
              CALL RANDOM_NUMBER(r)
              E = r * 10.0
              P_ads = 1.0 - SIN(theta) * (E / 10.0)
              CALL RANDOM_NUMBER(r)
              IF (P_ads > r) THEN
                  height(x,y) = height(x,y) + 1.0
                  P_diff = EXP(-ed / (kb * t_subs))
                  CALL RANDOM_NUMBER(r)
                  IF (P_diff > r) THEN
                      CALL RANDOM_NUMBER(r)
                      IF (r < 0.25 .AND. x > 1) x = x - 1
                      IF (r < 0.50 .AND. y > 1) y = y - 1
                      IF (r < 0.75 .AND. x < L) x = x + 1
                      IF (r < 1.00 .AND. y < L) y = y + 1
                      IF (x >= 1 .AND. x <= L .AND. y >= 1 .AND. y <= L) THEN
                          height(x,y) = height(x,y) + 1.0
                          height(x,y-1) = height(x,y-1) - 1.0
                      END IF
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

      SUBROUTINE save_data(height, L, filename)
          REAL, DIMENSION(L,L), INTENT(IN) :: height
          INTEGER, INTENT(IN) :: L
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER :: i, j
          OPEN(UNIT=10, FILE=filename, STATUS='REPLACE')
          DO i = 1, L
              DO j = 1, L
                  WRITE(10,*) i, j, height(i,j)
              END DO
          END DO
          CLOSE(10)
      END SUBROUTINE save_data
      END PROGRAM kmc_sputtering_simulation
