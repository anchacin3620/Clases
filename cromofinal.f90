PROGRAM simulacion_montecarlo_pulverizacion_3d_cr
    IMPLICIT NONE
    ! Parámetros generales
    INTEGER, PARAMETER :: TAM_REJILLA = 300, ALTURA_MAX = 150
    REAL, PARAMETER :: KB = 8.617E-5  ! Constante de Boltzmann (eV/K)
    REAL, PARAMETER :: TEMP_SUBSTRATO = 500.0  ! Temperatura del sustrato (K)
    REAL, PARAMETER :: TEMP_PLASMA = 1000.0    ! Temperatura del plasma (K)
    CHARACTER(LEN=2), PARAMETER :: MATERIAL_TYPE = 'Cr'  ! Material: Cromo
    REAL, PARAMETER :: ADHESION_BASE = 0.98     ! Probabilidad base de adhesión
    REAL, PARAMETER :: ENERGIA_DIFUSION_BASE = 0.6  ! Energía de difusión (eV)
    INTEGER, PARAMETER :: MAX_PART = 100000      ! Máximo número de partículas

    ! Arreglos
    INTEGER, DIMENSION(TAM_REJILLA, TAM_REJILLA, ALTURA_MAX) :: material
    INTEGER, DIMENSION(TAM_REJILLA, TAM_REJILLA) :: altura_max_local

    ! Variables
    INTEGER :: i, j, num_eventos, semilla(8), depositos, difusiones, num_corridas, corrida, modo
    REAL :: tiempo, prob_difusion_const, energia_difusion, tasa_eventos, adhesion, energia_promedio
    REAL :: rugosidad_rms
    CHARACTER(LEN=50) :: nombre_archivo

    ! Solicitar número de corridas
    PRINT *, 'Ingrese el número de corridas deseadas:'
    READ *, num_corridas
    IF (num_corridas <= 0) THEN
        PRINT *, 'Error: El número de corridas debe ser positivo'
        STOP
    END IF

    ! Bucle para múltiples corridas
    DO corrida = 1, num_corridas
        ! Bucle para diferentes modos de deposición
        DO modo = 1, 3
            ! Inicializar variables
            num_eventos = 0
            depositos = 0
            difusiones = 0
            tiempo = 0.0

            ! Configurar parámetros según el modo de deposición
            IF (modo == 1) THEN
                energia_promedio = 0.3  ! Evaporación térmica
                tasa_eventos = 10000.0
            ELSE IF (modo == 2) THEN
                energia_promedio = 3.5  ! Pulverización catódica
                tasa_eventos = 50000.0
            ELSE
                energia_promedio = 7.0  ! Láser pulsado
                tasa_eventos = 100000.0
            END IF

            ! Generar parámetros aleatorios
            CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(energia_difusion)
            energia_difusion = ENERGIA_DIFUSION_BASE - 0.1 + 0.2 * energia_difusion  ! Rango: ±0.1 eV
            CALL RANDOM_NUMBER(adhesion)
            adhesion = ADHESION_BASE - 0.03 + 0.06 * adhesion  ! Rango: ±0.03

            ! Precalcular constante de difusión
            prob_difusion_const = EXP(-energia_difusion / (KB * TEMP_SUBSTRATO))

            ! Inicializar semilla
            CALL SYSTEM_CLOCK(i)
            semilla = i + 37 * (corrida + modo)
            CALL RANDOM_SEED(PUT=semilla)

            ! Inicializar rejilla
            CALL inicializar_rejilla(material, altura_max_local, TAM_REJILLA, ALTURA_MAX)

            ! Simulación KMC
            DO WHILE (num_eventos < MAX_PART)
                CALL simular_evento(material, altura_max_local, num_eventos, MAX_PART, &
                                    tasa_eventos, tiempo, energia_difusion, TEMP_SUBSTRATO, &
                                    KB, adhesion, TEMP_PLASMA, TAM_REJILLA, ALTURA_MAX, &
                                    depositos, difusiones, prob_difusion_const, energia_promedio)
            END DO

            ! Calcular rugosidad
            CALL calcular_rugosidad(material, TAM_REJILLA, ALTURA_MAX, rugosidad_rms)

            ! Generar nombre de archivo único
            WRITE(nombre_archivo, '(A,I0.3,A,I0,A,F4.2,A,F6.0,A,F4.1,A)') &
                  MATERIAL_TYPE//'_3d_altura_c', corrida, '_m', modo, '_ed', energia_difusion, &
                  '_te', tasa_eventos, '_ep', energia_promedio, '.dat'

            ! Guardar resultados
            CALL guardar_datos(material, TAM_REJILLA, ALTURA_MAX, nombre_archivo, &
                               energia_difusion, TEMP_SUBSTRATO, tasa_eventos, adhesion, rugosidad_rms)

            PRINT *, 'Corrida ', corrida, ' Modo ', modo, ' Depósitos = ', depositos, &
                     ' Difusiones = ', difusiones, ' Tiempo = ', tiempo, ' Rugosidad RMS = ', rugosidad_rms
        END DO
    END DO

    PRINT *, 'Simulaciones para ', MATERIAL_TYPE, ' completadas'
END PROGRAM simulacion_montecarlo_pulverizacion_3d_cr

! Subrutina para inicializar la rejilla (Mg(001) como sustrato hcp)
SUBROUTINE inicializar_rejilla(material, altura_max_local, tam_rejilla, altura_max)
    INTEGER, INTENT(IN) :: tam_rejilla, altura_max
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla, altura_max), INTENT(OUT) :: material
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla), INTENT(OUT) :: altura_max_local
    INTEGER :: i, j

    material = 0
    DO i = 1, tam_rejilla
        DO j = 1, tam_rejilla
            material(i, j, 1) = 3  ! Sustrato Mg(001)
        END DO
    END DO
    altura_max_local = 1  ! Altura inicial del sustrato
END SUBROUTINE inicializar_rejilla

! Función para encontrar la altura máxima en una columna
INTEGER FUNCTION find_z_max(columna)
    INTEGER, INTENT(IN) :: columna(ALTURA_MAX)
    INTEGER :: k
    DO k = ALTURA_MAX, 1, -1
        IF (columna(k) /= 0) THEN
            find_z_max = k
            RETURN
        END IF
    END DO
    find_z_max = 0
END FUNCTION find_z_max

! Subrutina para generar números aleatorios seguros
SUBROUTINE safe_random(r)
    REAL, INTENT(OUT) :: r
    CALL RANDOM_NUMBER(r)
    r = 1E-10 + (1.0 - 2E-10) * r
END SUBROUTINE safe_random

! Subrutina para simular un evento (deposición o difusión)
SUBROUTINE simular_evento(material, altura_max_local, num_eventos, max_part, &
                          tasa_eventos, tiempo, energia_difusion, temp_substrato, &
                          kb, adhesion, temp_plasma, tam_rejilla, altura_max, &
                          depositos, difusiones, prob_difusion_const, energia_promedio)
    INTEGER, INTENT(IN) :: tam_rejilla, altura_max, max_part
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla, altura_max), INTENT(INOUT) :: material
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla), INTENT(INOUT) :: altura_max_local
    INTEGER, INTENT(INOUT) :: num_eventos, depositos, difusiones
    REAL, INTENT(IN) :: tasa_eventos, energia_difusion, temp_substrato, kb, adhesion, &
                        temp_plasma, prob_difusion_const, energia_promedio
    REAL, INTENT(INOUT) :: tiempo
    INTEGER :: x, y, z_max, nuevo_x, nuevo_y
    REAL :: energia, prob_adhesion, r, r_tiempo
    REAL, PARAMETER :: MIN_R = 1E-10, MAX_ENERGIA = 5.0
    INTEGER :: vecinos(6, 2)  ! Rejilla hexagonal
    INTEGER :: i, vecino_idx

    num_eventos = num_eventos + 1
    IF (num_eventos > max_part) RETURN

    ! Avanzar el tiempo
    CALL safe_random(r_tiempo)
    tiempo = tiempo + (-LOG(MAX(MIN_R, r_tiempo)) / tasa_eventos)

    ! Seleccionar sitio aleatorio
    CALL safe_random(r)
    x = 1 + FLOOR(r * REAL(tam_rejilla-1))
    CALL safe_random(r)
    y = 1 + FLOOR(r * REAL(tam_rejilla-1))
    z_max = MAXVAL(altura_max_local(MAX(1,x-1):MIN(tam_rejilla,x+1), &
                                    MAX(1,y-1):MIN(tam_rejilla,y+1)))

    ! Deposición
    CALL safe_random(r)
    energia = energia_promedio + (-0.5 + r) * MAX_ENERGIA
    prob_adhesion = adhesion * (1.0 + (energia / MAX_ENERGIA) * 0.2)
    IF (prob_adhesion > r .AND. z_max < altura_max) THEN
        depositos = depositos + 1
        material(x, y, z_max+1) = 2
        altura_max_local(x, y) = z_max + 1
    END IF

    ! Difusión (en rejilla hexagonal)
    CALL safe_random(r)
    IF (prob_difusion_const > r .AND. z_max >= 1) THEN
        ! Definir vecinos para rejilla hexagonal
        vecinos = RESHAPE([x+1,y, x-1,y, x,y+1, x,y-1, x+1,y-1, x-1,y+1], [6,2])
        CALL safe_random(r)
        vecino_idx = 1 + FLOOR(r * 6.0)
        nuevo_x = MODULO(vecinos(vecino_idx,1)-1, tam_rejilla) + 1
        nuevo_y = MODULO(vecinos(vecino_idx,2)-1, tam_rejilla) + 1
        IF (nuevo_x >= 1 .AND. nuevo_x <= tam_rejilla .AND. &
            nuevo_y >= 1 .AND. nuevo_y <= tam_rejilla) THEN
            IF (altura_max_local(nuevo_x, nuevo_y) < altura_max_local(x, y)) THEN
                difusiones = difusiones + 1
                material(nuevo_x, nuevo_y, altura_max_local(nuevo_x, nuevo_y)+1) = 2
                material(x, y, altura_max_local(x, y)) = 0
                altura_max_local(x, y) = find_z_max(material(x, y, :))
                altura_max_local(nuevo_x, nuevo_y) = find_z_max(material(nuevo_x, nuevo_y, :))
            END IF
        END IF
    END IF
END SUBROUTINE simular_evento

! Subrutina para calcular rugosidad RMS
SUBROUTINE calcular_rugosidad(material, tam_rejilla, altura_max, rms)
    INTEGER, INTENT(IN) :: tam_rejilla, altura_max
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla, altura_max), INTENT(IN) :: material
    REAL, INTENT(OUT) :: rms
    INTEGER :: i, j, z_max
    REAL :: altura_media, suma
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla) :: alturas

    ! Calcular alturas
    DO i = 1, tam_rejilla
        DO j = 1, tam_rejilla
            alturas(i, j) = find_z_max(material(i, j, :))
        END DO
    END DO

    ! Calcular altura media
    altura_media = SUM(alturas) / REAL(tam_rejilla * tam_rejilla)

    ! Calcular RMS
    suma = 0.0
    DO i = 1, tam_rejilla
        DO j = 1, tam_rejilla
            suma = suma + (REAL(alturas(i, j)) - altura_media)**2
        END DO
    END DO
    rms = SQRT(suma / REAL(tam_rejilla * tam_rejilla))
END SUBROUTINE calcular_rugosidad

! Subrutina para calcular densidad de núcleos (algoritmo Hoshen-Kopelman simplificado)
SUBROUTINE calcular_densidad_nucleos(material, tam_rejilla, altura_max, densidad)
    INTEGER, INTENT(IN) :: tam_rejilla, altura_max
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla, altura_max), INTENT(IN) :: material
    REAL, INTENT(OUT) :: densidad
    INTEGER :: i, j, z_max, num_nucleos
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla) :: etiquetas
    INTEGER :: etiqueta_actual, max_etiqueta
    INTEGER, ALLOCATABLE :: uniones(:)
    INTEGER :: x, y, vecino_idx, nuevo_x, nuevo_y
    INTEGER :: vecinos(6, 2)

    etiquetas = 0
    etiqueta_actual = 0
    ALLOCATE(uniones(tam_rejilla * tam_rejilla))

    ! Etiquetar núcleos en la capa superior
    DO i = 1, tam_rejilla
        DO j = 1, tam_rejilla
            z_max = find_z_max(material(i, j, :))
            IF (material(i, j, z_max) == 2) THEN
                ! Verificar vecinos (rejilla hexagonal)
                vecinos = RESHAPE([i+1,j, i-1,j, i,j+1, i,j-1, i+1,j-1, i-1,j+1], [6,2])
                max_etiqueta = 0
                DO vecino_idx = 1, 6
                    nuevo_x = MODULO(vecinos(vecino_idx,1)-1, tam_rejilla) + 1
                    nuevo_y = MODULO(vecinos(vecino_idx,2)-1, tam_rejilla) + 1
                    IF (nuevo_x >= 1 .AND. nuevo_x <= tam_rejilla .AND. &
                        nuevo_y >= 1 .AND. nuevo_y <= tam_rejilla) THEN
                        IF (etiquetas(nuevo_x, nuevo_y) > 0) THEN
                            max_etiqueta = MAX(max_etiqueta, etiquetas(nuevo_x, nuevo_y))
                        END IF
                    END IF
                END DO
                IF (max_etiqueta == 0) THEN
                    etiqueta_actual = etiqueta_actual + 1
                    etiquetas(i, j) = etiqueta_actual
                    uniones(etiqueta_actual) = etiqueta_actual
                ELSE
                    etiquetas(i, j) = max_etiqueta
                END IF
            END IF
        END DO
    END DO

    ! Contar núcleos únicos
    num_nucleos = 0
    DO i = 1, etiqueta_actual
        IF (uniones(i) == i) num_nucleos = num_nucleos + 1
    END DO
    densidad = REAL(num_nucleos) / REAL(tam_rejilla * tam_rejilla)
    DEALLOCATE(uniones)
END SUBROUTINE calcular_densidad_nucleos

! Subrutina para guardar datos
SUBROUTINE guardar_datos(material, tam_rejilla, altura_max, nombre_archivo, &
                         energia_difusion, temp_substrato, tasa_eventos, adhesion, rms)
    INTEGER, INTENT(IN) :: tam_rejilla, altura_max
    INTEGER, DIMENSION(tam_rejilla, tam_rejilla, altura_max), INTENT(IN) :: material
    CHARACTER(LEN=*), INTENT(IN) :: nombre_archivo
    REAL, INTENT(IN) :: energia_difusion, temp_substrato, tasa_eventos, adhesion, rms
    INTEGER :: i, j, z_max
    REAL :: densidad_nucleos

    CALL calcular_densidad_nucleos(material, tam_rejilla, altura_max, densidad_nucleos)

    OPEN(UNIT=10, FILE=nombre_archivo, STATUS='REPLACE')
    WRITE(10, '(A)') '# Simulación KMC 3D para pulverización de '//MATERIAL_TYPE
    WRITE(10, '(A, I0)') '# Tamaño de la rejilla: ', tam_rejilla
    WRITE(10, '(A, I0)') '# Altura máxima: ', altura_max
    WRITE(10, '(A, F8.2)') '# Energía de difusión (eV): ', energia_difusion
    WRITE(10, '(A, F8.1)') '# Temperatura del sustrato (K): ', temp_substrato
    WRITE(10, '(A, F10.0)') '# Tasa de eventos: ', tasa_eventos
    WRITE(10, '(A, F8.3)') '# Probabilidad de adhesión: ', adhesion
    WRITE(10, '(A, F8.3)') '# Rugosidad RMS: ', rms
    WRITE(10, '(A, F8.5)') '# Densidad de núcleos (1/área): ', densidad_nucleos
    WRITE(10, '(A)') '# Formato: i j altura(i,j)'
    DO i = 1, tam_rejilla
        DO j = 1, tam_rejilla
            z_max = find_z_max(material(i, j, :))
            WRITE(10, '(I4, I4, I4)') i, j, z_max
        END DO
    END DO
    CLOSE(10)
END SUBROUTINE guardar_datos
