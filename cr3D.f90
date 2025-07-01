      PROGRAM MonteCarloSputtering3D_Cr
      IMPLICIT NONE
      INTEGER, PARAMETER :: L = 300, H = 100, MAX_PART = 20000
      REAL, PARAMETER :: KB = 8.617E-5, ED = 0.6, T_SUBS = 500.0, T_PLASMA = 1000.0
      REAL, PARAMETER :: STICKING = 0.9, LAMBDA = 20000.0
      INTEGER, DIMENSION(L,L,H) :: material
      INTEGER :: i, j, k, n_events, seed
      REAL :: time

      CALL SYSTEM_CLOCK(seed)
      CALL RANDOM_SEED(seed)

      CALL initialize_grid(material, L, H)

      n_events = 0
      time = 0.0
      DO WHILE (n_events < MAX_PART)
          CALL simulate_event(material, n_events, MAX_PART, LAMBDA, time, ED, T_SUBS, KB, STICKING, T_PLASMA)
          IF (n_events >= MAX_PART) EXIT
      END DO

      CALL save_data(material, L, H, 'cr_3d_height.dat')

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

      SUBROUTINE simulate_event(material, n_events, max_part, lambda, time, ed, t_subs, kb, sticking, t_plasma)
          INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: material
          INTEGER, INTENT(INOUT) :: n_events
          INTEGER, INTENT(IN) :: max_part
          REAL, INTENT(IN) :: lambda, time, ed, t_subs, kb, sticking, t_plasma
          INTEGER :: x, y, z_max, i, j
          REAL :: E, P_stick, r, P_diff, r_time
          REAL, PARAMETER :: pi = 3.14159265359

          CALL RANDOM_NUMBER(r)
          n_events = n_events + FLOOR(-LOG(r) * lambda)
          IF (n_events > max_part) n_events = max_part

          CALL RANDOM_NUMBER(r_time)
          time = time - LOG(r_time) / lambda

          IF (n_events <= max_part) THEN
              CALL RANDOM_NUMBER(r)
              x = MIN(MAX(1, INT(r * L) + 1), L)
              CALL RANDOM_NUMBER(r)
              y = MIN(MAX(1, INT(r * L) + 1), L)
              z_max = 2
              DO i = MAX(1,x-1), MIN(L,x+1)
                  DO j = MAX(1,y-1), MIN(L,y+1)
                      z_max = MAX(z_max, MAXVAL(material(i,j,1:H)))
                  END DO
              END DO
              CALL RANDOM_NUMBER(r)
              E = -KB * t_plasma * LOG(r)
              P_stick = sticking * EXP(-E / (KB * t_plasma))
              CALL RANDOM_NUMBER(r)
              IF (P_stick > r) THEN
                  material(x,y,z_max+1) = 1
                  P_diff = EXP(-ed / (kb * t_subs))
                  CALL RANDOM_NUMBER(r)
                  IF (P_diff > r .AND. z_max > 2) THEN
                      CALL RANDOM_NUMBER(r)
                      IF (r < 0.25 .AND. x > 1) x = x - 1
                      IF (r < 0.50 .AND. y > 1) y = y - 1
                      IF (r < 0.75 .AND. x < L) x = x + 1
                      IF (r < 1.00 .AND. y < L) y = y + 1
                      IF (x >= 1 .AND. x <= L .AND. y >= 1 .AND. y <= L) THEN
                          material(x,y,z_max+1) = 1
                          material(x,y-1,z_max+1) = 0
                      END IF
                  END IF
              END IF
          END IF
      END SUBROUTINE simulate_event

      SUBROUTINE save_data(material, L, H, filename)
          INTEGER, DIMENSION(L,L,H), INTENT(IN) :: material
          INTEGER, INTENT(IN) :: L, H
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER :: i, j, z_max
          OPEN(UNIT=10, FILE=filename, STATUS='REPLACE')
          DO i = 1, L
              DO j = 1, L
                  z_max = MAXVAL(material(i,j,1:H))
                  WRITE(10,*) i, j, z_max
              END DO
          END DO
          CLOSE(10)
      END SUBROUTINE save_data
      END PROGRAM MonteCarloSputtering3D_Cr
