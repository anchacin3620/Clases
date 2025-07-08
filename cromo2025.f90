! Programa de Simulación Monte Carlo elaborado por Ángel Chacín Ávila, codificado en GNU Fortran y bajo licencia GPL
PROGRAM simulacion_montecarlo_pulverizacion_3d_cr
    IMPLICIT NONE
    INTEGER, PARAMETER :: TAM_REJILLA = 300, ALTURA_MAX = 100, MAX_PARTICULAS = 10000
    REAL, PARAMETER :: KB = 8.617E-5, ENERGIA_DIFUSION = 0.6, TEMP_SUBSTRATO = 500.0, TEMP_PLASMA = 1000.0
    REAL, PARAMETER :: ADHESION = 0.9, TASA_EVENTOS = 3000000.0
    INTEGER, DIMENSION(TAM_REJILLA,TAM_REJILLA,ALTURA_MAX) :: material
    INTEGER, DIMENSION(TAM_REJILLA,TAM_REJILLA) :: altura_max_local
    INTEGER :: i, j, k, num_eventos, semilla(8), depositos
    REAL :: tiempo, prob_difusion_const

    ! Validar parámetros físicos
    IF (ENERGIA_DIFUSION <= 0.0 .OR. TEMP_SUBSTRATO <= 0.0 .OR. TEMP_PLASMA <= 0.0 .OR. TASA_EVENTOS <= 0.0 .OR. ADHESION <= 0.0) THEN
        PRINT *, 'Error: ENERGIA_DIFUSION, TEMP_SUBSTRATO, TEMP_PLASMA, TASA_EVENTOS y ADHESION deben ser positivos'
        STOP
    END IF

    ! Precalcular constante de difusión
    prob_difusion_const = EXP(-ENERGIA_DIFUSION / (KB * TEMP_SUBSTRATO))

    ! Inicializar semilla para números aleatorios
    CALL SYSTEM_CLOCK(i)
    semilla = i + 37 * (/ (j - 1, j = 1, 8) /)
    CALL RANDOM_SEED(PUT=semilla)

    ! Inicializar rejilla
    CALL inicializar_rejilla(material, altura_max_local, TAM_REJILLA, ALTURA_MAX)

    ! Simulación con distribución de Poisson
    num_eventos = 0
    depositos = 0
    tiempo = 0.0
    DO WHILE (num_eventos < MAX_PARTICULAS)
        CALL simular_evento(material, altura_max_local, num_eventos, MAX_PARTICULAS, TASA_EVENTOS, tiempo, ENERGIA_DIFUSION, TEMP_SUBSTRATO, KB, ADHESION, TEMP_PLASMA, TAM_REJILLA, ALTURA_MAX, depositos, prob_difusion_const)
        IF (num_eventos >= MAX_PARTICULAS) EXIT
    END DO
    PRINT *, 'Total de depósitos: ', depositos

    ! Guardar resultados para Gnuplot
    CALL guardar_datos(material, TAM_REJILLA, ALTURA_MAX, 'cr_3d_altura.dat', ENERGIA_DIFUSION, TEMP_SUBSTRATO, TASA_EVENTOS, ADHESION)

    PRINT *, 'Simulación 3D para Cr completada. Ver cr_3d_altura.dat'
    STOP

CONTAINS

    SUBROUTINE inicializar_rejilla(material, altura_max_local, tam_rejilla, altura_max)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(OUT) :: material
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla), INTENT(OUT) :: altura_max_local
        INTEGER, INTENT(IN) :: tam_rejilla, altura_max
        INTEGER :: i, j
        material = 0
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                material(i,j,1) = 3
                material(i,j,2) = 3
            END DO
        END DO
        altura_max_local = 2  ! Capa base del sustrato
    END SUBROUTINE inicializar_rejilla

    INTEGER FUNCTION find_z_max(columna)
        INTEGER, INTENT(IN) :: columna(ALTURA_MAX)
        INTEGER :: k
        DO k = ALTURA_MAX, 1, -1
            IF (columna(k) > 0) THEN
                find_z_max = k
                RETURN
            END IF
        END DO
        find_z_max = 0
    END FUNCTION find_z_max

    SUBROUTINE safe_random(r)
        REAL, INTENT(OUT) :: r
        CALL RANDOM_NUMBER(r)
        r = 1E-10 + (1 - 2E-10) * r
    END SUBROUTINE safe_random

    SUBROUTINE simular_evento(material, altura_max_local, num_eventos, max_part, tasa_eventos, tiempo, energia_difusion, temp_substrato, kb, adhesion, temp_plasma, tam_rejilla, altura_max, depositos, prob_difusion_const)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(INOUT) :: material
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla), INTENT(INOUT) :: altura_max_local
        INTEGER, INTENT(INOUT) :: num_eventos, depositos
        INTEGER, INTENT(IN) :: max_part, tam_rejilla, altura_max
        REAL, INTENT(IN) :: tasa_eventos, energia_difusion, temp_substrato, kb, adhesion, temp_plasma, prob_difusion_const
        REAL, INTENT(INOUT) :: tiempo
        INTEGER :: x, y, z_max, nuevo_x, nuevo_y, direction
        REAL :: energia, prob_adhesion, r, r_tiempo
        REAL, PARAMETER :: MIN_R = 1.0E-10, MAX_ENERGIA = 5.0

        CALL RANDOM_NUMBER(r)
        num_eventos = num_eventos + INT(MAX(1.0, FLOOR(-LOG(MAX(MIN_R, r)) * tasa_eventos)))  ! Corrección en línea 92
        IF (num_eventos > max_part) num_eventos = max_part

        CALL RANDOM_NUMBER(r_tiempo)
        r_tiempo = MAX(MIN_R, r_tiempo)
        tiempo = tiempo + (-LOG(r_tiempo) / tasa_eventos)

        IF (num_eventos <= max_part) THEN
            CALL safe_random(r)
            x = 1 + FLOOR(r * (tam_rejilla-1))
            CALL safe_random(r)
            y = 1 + FLOOR(r * (tam_rejilla-1))
            z_max = MAXVAL(altura_max_local(MAX(1,x-1):MIN(tam_rejilla,x+1), MAX(1,y-1):MIN(tam_rejilla,y+1)))
            CALL safe_random(r)
            energia = r * MAX_ENERGIA
            CALL safe_random(r)
            prob_adhesion = adhesion * (1.0 - r * 0.5)
            ! PRINT *, 'Evento en (', x, ',', y, '), z_max = ', z_max, ', prob_adhesion = ', prob_adhesion, ', r = ', r  ! Comentar en producción
            IF (prob_adhesion > r .AND. z_max < altura_max) THEN
                depositos = depositos + 1
                material(x,y,z_max+1) = 1
                altura_max_local(x,y) = z_max + 1
                ! PRINT *, 'Depósito en (', x, ',', y, '), nuevo z_max = ', z_max + 1  ! Comentar en producción
                IF (z_max > 2) THEN
                    CALL safe_random(r)
                    IF (prob_difusion_const > r) THEN
                        nuevo_x = MODULO(x + MERGE(-1,1,r<0.5) - 1, tam_rejilla) + 1
                        nuevo_y = MODULO(y + MERGE(-1,1,r>0.5) - 1, tam_rejilla) + 1
                        material(nuevo_x, nuevo_y, z_max) = 1
                        material(x, y, z_max) = 0
                        altura_max_local(x,y) = find_z_max(material(x,y,:))
                        altura_max_local(nuevo_x,nuevo_y) = find_z_max(material(nuevo_x,nuevo_y,:))
                    END IF
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
        WRITE(10, '(A)') '# Simulación KMC 3D para pulverización (Cr)'
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
                z_max = find_z_max(material(i,j,:))
                WRITE(10, '(I4, I4, I4)') i, j, z_max
            END DO
        END DO
        CLOSE(10)
    END SUBROUTINE guardar_datos
END PROGRAM simulacion_montecarlo_pulverizacion_3d_cr
