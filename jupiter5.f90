PROGRAM simulation_kmc_unificado
    IMPLICIT NONE
    ! Parámetros constantes (reducidos para pruebas)
    INTEGER, PARAMETER :: TAM_REJILLA = 50, ALTURA_MAX = 30
    REAL, PARAMETER :: KB = 8.617E-5, TEMP_SUBSTRATO = 500.0, TEMP_PLASMA = 1000.0
    
    ! Parámetros específicos por material
    REAL :: ADHESION_BASE, ENERGIA_DIFUSION_MIN, ENERGIA_DIFUSION_MAX
    INTEGER :: MATERIAL_ID
    
    ! Variables principales
    INTEGER, DIMENSION(TAM_REJILLA,TAM_REJILLA,ALTURA_MAX) :: material
    INTEGER, DIMENSION(TAM_REJILLA,TAM_REJILLA) :: altura_max_local
    INTEGER :: i, j, k, num_eventos, semilla(8), depositos, difusiones, num_corridas, corrida
    REAL :: tiempo, prob_difusion_const, energia_difusion, tasa_eventos, adhesion, energia_promedio
    CHARACTER(LEN=50) :: nombre_archivo, material_name
    INTEGER :: num_snapshots, snapshot_actual, input_value
    INTEGER, PARAMETER :: INTERVALO_SNAPSHOT = 1000, MAX_EVENTOS = 5000  ! Reducido para pruebas
    
    ! Variables para métricas avanzadas
    REAL :: densidad_islas, tamano_promedio_islas, cobertura_fraccional, rugosidad_rms
    REAL :: anisotropia_superficial, longitud_difusion_promedio
    INTEGER :: num_islas, total_atomos_superficie
    
    ! Variables para modos de deposición
    INTEGER :: modo_deposicion
    REAL :: energia_thermal = 2.0, energia_ballistic = 5.0
    REAL :: tasa_thermal = 50000.0, tasa_ballistic = 20000.0
    
    ! Selección de material
    PRINT *, 'Seleccione el material (1=Cr, 2=Fe):'
    READ *, MATERIAL_ID
    IF (MATERIAL_ID == 1) THEN
        material_name = 'Cr'
        ADHESION_BASE = 0.98
        ENERGIA_DIFUSION_MIN = 0.01
        ENERGIA_DIFUSION_MAX = 0.25
    ELSE IF (MATERIAL_ID == 2) THEN
        material_name = 'Fe'
        ADHESION_BASE = 0.95
        ENERGIA_DIFUSION_MIN = 0.5
        ENERGIA_DIFUSION_MAX = 1.0
    ELSE
        PRINT *, 'Material no válido'
        STOP
    END IF
    
    ! Selección de modo de deposición
    PRINT *, 'Modo de deposición (1=Térmico, 2=Balístico, 3=Aleatorio):'
    READ *, modo_deposicion
    
    ! Solicitar número de corridas
    PRINT *, 'Ingrese el número de corridas deseadas:'
    READ *, num_corridas
    IF (num_corridas <= 0) THEN
        PRINT *, 'Error: El número de corridas debe ser positivo'
        STOP
    END IF
    
    ! Inicializar snapshots
    PRINT *, '¿Es la primera corrida? (1=Sí, 0=No):'
    READ *, input_value
    IF (input_value == 1) THEN
        num_snapshots = 0
    ELSE
        PRINT *, 'Ingrese el número del último snapshot:'
        READ *, num_snapshots
    END IF
    
    ! Bucle para múltiples corridas
    DO corrida = 1, num_corridas
        ! Generar parámetros aleatorios según modo de deposición
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(energia_difusion)
        energia_difusion = ENERGIA_DIFUSION_MIN + (ENERGIA_DIFUSION_MAX - ENERGIA_DIFUSION_MIN) * energia_difusion
        
        ! Configurar parámetros según modo de deposición
        IF (modo_deposicion == 1) THEN  ! Térmico
            energia_promedio = energia_thermal
            tasa_eventos = tasa_thermal
        ELSE IF (modo_deposicion == 2) THEN  ! Balístico
            energia_promedio = energia_ballistic
            tasa_eventos = tasa_ballistic
        ELSE  ! Aleatorio
            CALL RANDOM_NUMBER(tasa_eventos)
            tasa_eventos = 10000.0 + 90000.0 * tasa_eventos
            CALL RANDOM_NUMBER(energia_promedio)
            energia_promedio = 2.0 + 3.0 * energia_promedio
        END IF
        
        CALL RANDOM_NUMBER(adhesion)
        adhesion = ADHESION_BASE + 0.03 * adhesion
        
        ! Precalcular constante de difusión
        prob_difusion_const = EXP(-energia_difusion / (KB * TEMP_SUBSTRATO))
        
        ! Inicializar semilla para números aleatorios
        CALL SYSTEM_CLOCK(i)
        semilla = i + 37 * (/ (j - 1, j = 1, 8) /) + corrida
        CALL RANDOM_SEED(PUT=semilla)
        
        ! Inicializar rejilla
        CALL inicializar_rejilla(material, altura_max_local, TAM_REJILLA, ALTURA_MAX)
        
        ! Simulación
        num_eventos = 0
        depositos = 0
        difusiones = 0
        tiempo = 0.0
        snapshot_actual = 0
        
        DO WHILE (num_eventos < MAX_EVENTOS)
            CALL simular_evento(material, altura_max_local, num_eventos, MAX_EVENTOS, &
                tasa_eventos, tiempo, energia_difusion, TEMP_SUBSTRATO, &
                KB, adhesion, TEMP_PLASMA, TAM_REJILLA, ALTURA_MAX, &
                depositos, difusiones, prob_difusion_const, energia_promedio, MATERIAL_ID)
                
            ! Guardar snapshot periódicamente
            IF (MOD(num_eventos, INTERVALO_SNAPSHOT) == 0 .AND. num_eventos > 0) THEN
                snapshot_actual = snapshot_actual + 1
                num_snapshots = num_snapshots + 1
                
                ! Calcular métricas avanzadas
                CALL calcular_mecanismos(material, altura_max_local, TAM_REJILLA, ALTURA_MAX, &
                    densidad_islas, tamano_promedio_islas, num_islas, total_atomos_superficie)
                    
                CALL calcular_propiedades_fisicas(material, altura_max_local, TAM_REJILLA, ALTURA_MAX, &
                    cobertura_fraccional, rugosidad_rms, anisotropia_superficial)
                    
                CALL calcular_longitud_difusion(difusiones, tiempo, longitud_difusion_promedio)
                
                WRITE(nombre_archivo, '(A,A1,A,I0.3,A,I0.3,A,F4.2,A,F6.0,A,F4.1,A)') &
                    TRIM(material_name), '_', '3d_altura_c', corrida, '_s', snapshot_actual, '_ed', energia_difusion, &
                    '_te', tasa_eventos, '_ep', energia_promedio, '.dat'
                    
                CALL guardar_datos(material, TAM_REJILLA, ALTURA_MAX, nombre_archivo, &
                    energia_difusion, TEMP_SUBSTRATO, tasa_eventos, adhesion, energia_promedio, num_eventos, &
                    densidad_islas, tamano_promedio_islas, cobertura_fraccional, rugosidad_rms, &
                    anisotropia_superficial, longitud_difusion_promedio, num_islas, total_atomos_superficie)
            END IF
            
            IF (num_eventos >= MAX_EVENTOS) EXIT
        END DO
        
        ! Guardar snapshot final
        snapshot_actual = snapshot_actual + 1
        num_snapshots = num_snapshots + 1
        
        ! Calcular métricas finales
        CALL calcular_mecanismos(material, altura_max_local, TAM_REJILLA, ALTURA_MAX, &
            densidad_islas, tamano_promedio_islas, num_islas, total_atomos_superficie)
            
        CALL calcular_propiedades_fisicas(material, altura_max_local, TAM_REJILLA, ALTURA_MAX, &
            cobertura_fraccional, rugosidad_rms, anisotropia_superficial)
            
        CALL calcular_longitud_difusion(difusiones, tiempo, longitud_difusion_promedio)
        
        WRITE(nombre_archivo, '(A,A1,A,I0.3,A,I0.3,A,F4.2,A,F6.0,A,F4.1,A)') &
            TRIM(material_name), '_', '3d_altura_c', corrida, '_s', snapshot_actual, '_ed', energia_difusion, &
            '_te', tasa_eventos, '_ep', energia_promedio, '.dat'
            
        CALL guardar_datos(material, TAM_REJILLA, ALTURA_MAX, nombre_archivo, &
            energia_difusion, TEMP_SUBSTRATO, tasa_eventos, adhesion, energia_promedio, num_eventos, &
            densidad_islas, tamano_promedio_islas, cobertura_fraccional, rugosidad_rms, &
            anisotropia_superficial, longitud_difusion_promedio, num_islas, total_atomos_superficie)
            
        PRINT *, 'Corrida', corrida, ': Depósitos =', depositos, 'Difusiones =', difusiones, 'Tiempo =', tiempo
        PRINT *, 'Densidad de islas:', densidad_islas, 'Tamaño promedio:', tamano_promedio_islas
        PRINT *, 'Cobertura:', cobertura_fraccional, 'Rugosidad:', rugosidad_rms
        PRINT *, 'Anisotropía:', anisotropia_superficial, 'Longitud difusión:', longitud_difusion_promedio
    END DO
    
    PRINT *, 'Simulaciones para', TRIM(material_name), 'completadas. Total de snapshots:', num_snapshots
    PRINT *, 'Elaborado por Ángel Chacín Ávila, se publica con licencia GPL'
    STOP

CONTAINS

    ! Inicializar rejilla con sustrato
    SUBROUTINE inicializar_rejilla(material, altura_max_local, tam_rejilla, altura_max)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(OUT) :: material
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla), INTENT(OUT) :: altura_max_local
        INTEGER, INTENT(IN) :: tam_rejilla, altura_max
        INTEGER :: i, j
        material = 0
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                material(i,j,1) = -1  ! Sustrato (valor negativo)
                material(i,j,2) = -1  ! Sustrato (valor negativo)
            END DO
        END DO
        altura_max_local = 2
    END SUBROUTINE inicializar_rejilla

    ! Encontrar altura máxima en una columna
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

    ! Generador de números aleatorios seguro
    SUBROUTINE safe_random(r)
        REAL, INTENT(OUT) :: r
        CALL RANDOM_NUMBER(r)
        r = 1E-10 + (1 - 2E-10) * r
    END SUBROUTINE safe_random

    ! Simular evento (depósito o difusión)
    SUBROUTINE simular_evento(material, altura_max_local, num_eventos, max_part, &
            tasa_eventos, tiempo, energia_difusion, temp_substrato, &
            kb, adhesion, temp_plasma, tam_rejilla, altura_max, &
            depositos, difusiones, prob_difusion_const, energia_promedio, material_id)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(INOUT) :: material
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla), INTENT(INOUT) :: altura_max_local
        INTEGER, INTENT(INOUT) :: num_eventos, depositos, difusiones
        INTEGER, INTENT(IN) :: max_part, tam_rejilla, altura_max, material_id
        REAL, INTENT(IN) :: tasa_eventos, energia_difusion, temp_substrato, kb, adhesion, temp_plasma, prob_difusion_const, energia_promedio
        REAL, INTENT(INOUT) :: tiempo
        INTEGER :: x, y, z_max, nuevo_x, nuevo_y
        REAL :: energia, prob_adhesion, r, r_tiempo, total_rate
        REAL, PARAMETER :: MIN_R = 1E-10, MAX_ENERGIA = 5.0
        
        num_eventos = num_eventos + 1
        IF (num_eventos > max_part) num_eventos = max_part
        
        ! Avanzar tiempo
        CALL RANDOM_NUMBER(r_tiempo)
        r_tiempo = MAX(MIN_R, r_tiempo)
        tiempo = tiempo + (-LOG(r_tiempo) / tasa_eventos)
        
        IF (num_eventos <= max_part) THEN
            ! Calcular tasas de eventos
            total_rate = tasa_eventos + prob_difusion_const * (tam_rejilla * tam_rejilla)
            CALL safe_random(r)
            
            IF (r < tasa_eventos / total_rate) THEN
                ! Evento de depósito
                CALL safe_random(r)
                x = 1 + FLOOR(r * (tam_rejilla-1))
                CALL safe_random(r)
                y = 1 + FLOOR(r * (tam_rejilla-1))
                z_max = MAXVAL(altura_max_local(MAX(1,x-1):MIN(tam_rejilla,x+1), &
                                MAX(1,y-1):MIN(tam_rejilla,y+1)))
                
                CALL safe_random(r)
                energia = energia_promedio + (r - 0.5) * MAX_ENERGIA
                CALL safe_random(r)
                prob_adhesion = adhesion * (1.0 + (energia / MAX_ENERGIA) * 0.2)
                
                IF (prob_adhesion > r .AND. z_max < altura_max) THEN
                    depositos = depositos + 1
                    material(x,y,z_max+1) = material_id  ! 1=Cr, 2=Fe
                    altura_max_local(x,y) = z_max + 1
                END IF
            ELSE
                ! Evento de difusión general
                CALL safe_random(r)
                x = 1 + FLOOR(r * (tam_rejilla-1))
                CALL safe_random(r)
                y = 1 + FLOOR(r * (tam_rejilla-1))
                z_max = altura_max_local(x,y)
                
                ! Solo difunden átomos depositados (no sustrato)
                IF (z_max > 2 .AND. material(x,y,z_max) > 0) THEN
                    CALL safe_random(r)
                    IF (prob_difusion_const > r) THEN
                        ! Elegir dirección de difusión
                        IF (r < 0.25) THEN
                            nuevo_x = MODULO(x - 2 + tam_rejilla - 1, tam_rejilla) + 1
                            nuevo_y = y
                        ELSE IF (r < 0.5) THEN
                            nuevo_x = MODULO(x, tam_rejilla) + 1
                            nuevo_y = y
                        ELSE IF (r < 0.75) THEN
                            nuevo_x = x
                            nuevo_y = MODULO(y - 2 + tam_rejilla - 1, tam_rejilla) + 1
                        ELSE
                            nuevo_x = x
                            nuevo_y = MODULO(y, tam_rejilla) + 1
                        END IF
                        
                        ! Verificar si el movimiento es válido
                        IF (nuevo_x >= 1 .AND. nuevo_x <= tam_rejilla .AND. &
                            nuevo_y >= 1 .AND. nuevo_y <= tam_rejilla) THEN
                            IF (material(nuevo_x, nuevo_y, z_max) == 0 .AND. &
                                altura_max_local(nuevo_x, nuevo_y) == z_max - 1) THEN
                                ! Realizar difusión
                                material(nuevo_x, nuevo_y, z_max) = material(x,y,z_max)
                                material(x,y,z_max) = 0
                                altura_max_local(x,y) = z_max - 1
                                altura_max_local(nuevo_x, nuevo_y) = z_max
                                difusiones = difusiones + 1
                            END IF
                        END IF
                    END IF
                END IF
            END IF
        END IF
    END SUBROUTINE simular_evento

    ! Calcular mecanismos de deposición (densidad y tamaño de islas)
    SUBROUTINE calcular_mecanismos(material, altura_max_local, tam_rejilla, altura_max, &
            densidad_islas, tamano_promedio_islas, num_islas, total_atomos_superficie)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(IN) :: material
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla), INTENT(IN) :: altura_max_local
        INTEGER, INTENT(IN) :: tam_rejilla, altura_max
        REAL, INTENT(OUT) :: densidad_islas, tamano_promedio_islas
        INTEGER, INTENT(OUT) :: num_islas, total_atomos_superficie
        
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla) :: visitado, altura_superficie
        INTEGER :: i, j, k, ii, jj, kk, tamano_isla
        INTEGER, DIMENSION(500000, 2) :: cola  ! Tamaño aumentado significativamente
        INTEGER :: frente, final
        INTEGER :: vecinos(4,2)
        INTEGER :: nuevo_x, nuevo_y
        
        ! Definir vecinos (arriba, abajo, izquierda, derecha)
        vecinos(1,1) = -1; vecinos(1,2) = 0
        vecinos(2,1) = 1; vecinos(2,2) = 0
        vecinos(3,1) = 0; vecinos(3,2) = -1
        vecinos(4,1) = 0; vecinos(4,2) = 1
        
        ! Inicializar
        visitado = 0
        altura_superficie = 0
        num_islas = 0
        total_atomos_superficie = 0
        
        ! Obtener altura de superficie
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                altura_superficie(i,j) = find_z_max(material(i,j,:))
            END DO
        END DO
        
        ! Encontrar islas usando BFS
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                IF (visitado(i,j) == 0 .AND. material(i,j,altura_superficie(i,j)) > 0) THEN
                    num_islas = num_islas + 1
                    tamano_isla = 0
                    
                    ! Inicializar BFS
                    cola(1,1) = i
                    cola(1,2) = j
                    frente = 1
                    final = 1
                    visitado(i,j) = 1
                    
                    DO WHILE (frente <= final .AND. frente <= 500000)
                        ii = cola(frente,1)
                        jj = cola(frente,2)
                        frente = frente + 1
                        tamano_isla = tamano_isla + 1
                        
                        ! Explorar vecinos (4-conectividad)
                        DO kk = 1, 4
                            nuevo_x = ii + vecinos(kk,1)
                            nuevo_y = jj + vecinos(kk,2)
                            
                            ! Aplicar condiciones de contorno periódicas
                            IF (nuevo_x < 1) nuevo_x = tam_rejilla
                            IF (nuevo_x > tam_rejilla) nuevo_x = 1
                            IF (nuevo_y < 1) nuevo_y = tam_rejilla
                            IF (nuevo_y > tam_rejilla) nuevo_y = 1
                            
                            ! Verificación adicional de límites
                            IF (nuevo_x >= 1 .AND. nuevo_x <= tam_rejilla .AND. &
                                nuevo_y >= 1 .AND. nuevo_y <= tam_rejilla) THEN
                                
                                IF (visitado(nuevo_x,nuevo_y) == 0 .AND. &
                                    material(nuevo_x,nuevo_y,altura_superficie(nuevo_x,nuevo_y)) > 0) THEN
                                    final = final + 1
                                    IF (final > 500000) THEN
                                        PRINT *, 'Advertencia: Cola llena en BFS'
                                        EXIT
                                    END IF
                                    cola(final,1) = nuevo_x
                                    cola(final,2) = nuevo_y
                                    visitado(nuevo_x,nuevo_y) = 1
                                END IF
                            END IF
                        END DO
                        IF (final > 500000) EXIT
                    END DO
                    
                    total_atomos_superficie = total_atomos_superficie + tamano_isla
                END IF
            END DO
        END DO
        
        ! Calcular métricas
        densidad_islas = REAL(num_islas) / (tam_rejilla * tam_rejilla)
        IF (num_islas > 0) THEN
            tamano_promedio_islas = REAL(total_atomos_superficie) / num_islas
        ELSE
            tamano_promedio_islas = 0.0
        END IF
    END SUBROUTINE calcular_mecanismos

    ! Calcular propiedades físicas
    SUBROUTINE calcular_propiedades_fisicas(material, altura_max_local, tam_rejilla, altura_max, &
            cobertura_fraccional, rugosidad_rms, anisotropia_superficial)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(IN) :: material
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla), INTENT(IN) :: altura_max_local
        INTEGER, INTENT(IN) :: tam_rejilla, altura_max
        REAL, INTENT(OUT) :: cobertura_fraccional, rugosidad_rms, anisotropia_superficial
        
        INTEGER :: i, j
        REAL :: altura_media, diff_x, diff_y
        INTEGER, DIMENSION(tam_rejilla, tam_rejilla) :: altura_superficie
        
        ! Obtener altura de superficie
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                altura_superficie(i,j) = find_z_max(material(i,j,:))
            END DO
        END DO
        
        ! Calcular altura media
        altura_media = SUM(REAL(altura_superficie)) / (tam_rejilla * tam_rejilla)
        
        ! Calcular rugosidad RMS
        rugosidad_rms = 0.0
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                rugosidad_rms = rugosidad_rms + (altura_superficie(i,j) - altura_media)**2
            END DO
        END DO
        rugosidad_rms = SQRT(rugosidad_rms / (tam_rejilla * tam_rejilla))
        
        ! Calcular cobertura fraccional
        cobertura_fraccional = 0.0
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                IF (altura_superficie(i,j) > 2) THEN
                    cobertura_fraccional = cobertura_fraccional + 1.0
                END IF
            END DO
        END DO
        cobertura_fraccional = cobertura_fraccional / (tam_rejilla * tam_rejilla)
        
        ! Calcular anisotropía superficial
        diff_x = 0.0
        diff_y = 0.0
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                diff_x = diff_x + ABS(altura_superficie(i,j) - altura_superficie(MODULO(i, tam_rejilla) + 1, j))
                diff_y = diff_y + ABS(altura_superficie(i,j) - altura_superficie(i, MODULO(j, tam_rejilla) + 1))
            END DO
        END DO
        diff_x = diff_x / (tam_rejilla * tam_rejilla)
        diff_y = diff_y / (tam_rejilla * tam_rejilla)
        
        IF (diff_y > 0.0) THEN
            anisotropia_superficial = diff_x / diff_y
        ELSE
            anisotropia_superficial = 1.0
        END IF
    END SUBROUTINE calcular_propiedades_fisicas

    ! Calcular longitud de difusión promedio
    SUBROUTINE calcular_longitud_difusion(difusiones, tiempo, longitud_difusion_promedio)
        INTEGER, INTENT(IN) :: difusiones
        REAL, INTENT(IN) :: tiempo
        REAL, INTENT(OUT) :: longitud_difusion_promedio
        
        IF (difusiones > 0 .AND. tiempo > 0.0) THEN
            longitud_difusion_promedio = SQRT(2.0 * REAL(difusiones) / tiempo)
        ELSE
            longitud_difusion_promedio = 0.0
        END IF
    END SUBROUTINE calcular_longitud_difusion

    ! Guardar datos con métricas avanzadas
    SUBROUTINE guardar_datos(material, tam_rejilla, altura_max, nombre_archivo, &
            energia_difusion, temp_substrato, tasa_eventos, adhesion, energia_promedio, num_eventos, &
            densidad_islas, tamano_promedio_islas, cobertura_fraccional, rugosidad_rms, &
            anisotropia_superficial, longitud_difusion_promedio, num_islas, total_atomos_superficie)
        INTEGER, DIMENSION(tam_rejilla,tam_rejilla,altura_max), INTENT(IN) :: material
        INTEGER, INTENT(IN) :: tam_rejilla, altura_max, num_eventos, num_islas, total_atomos_superficie
        CHARACTER(LEN=*), INTENT(IN) :: nombre_archivo
        REAL, INTENT(IN) :: energia_difusion, temp_substrato, tasa_eventos, adhesion, energia_promedio
        REAL, INTENT(IN) :: densidad_islas, tamano_promedio_islas, cobertura_fraccional
        REAL, INTENT(IN) :: rugosidad_rms, anisotropia_superficial, longitud_difusion_promedio
        
        INTEGER :: i, j, z_max
        
        OPEN(UNIT=10, FILE=nombre_archivo, STATUS='REPLACE')
        WRITE(10, '(A)') '# Simulación KMC 3D unificada'
        WRITE(10, '(A, I0)') '# Tamaño de la rejilla: ', tam_rejilla
        WRITE(10, '(A, I0)') '# Altura máxima: ', altura_max
        WRITE(10, '(A, I0)') '# Número de eventos: ', num_eventos
        WRITE(10, '(A, F8.3)') '# Energía de difusión (eV): ', energia_difusion
        WRITE(10, '(A, F8.1)') '# Temperatura del sustrato (K): ', temp_substrato
        WRITE(10, '(A, F8.1)') '# Tasa de eventos: ', tasa_eventos
        WRITE(10, '(A, F8.3)') '# Probabilidad de adhesión: ', adhesion
        WRITE(10, '(A, F8.1)') '# Energía promedio (eV): ', energia_promedio
        WRITE(10, '(A, F8.3)') '# Densidad de islas: ', densidad_islas
        WRITE(10, '(A, F8.1)') '# Tamaño promedio de islas: ', tamano_promedio_islas
        WRITE(10, '(A, F8.3)') '# Cobertura fraccional: ', cobertura_fraccional
        WRITE(10, '(A, F8.3)') '# Rugosidad RMS (capas): ', rugosidad_rms
        WRITE(10, '(A, F8.3)') '# Anisotropía superficial: ', anisotropia_superficial
        WRITE(10, '(A, F8.3)') '# Longitud de difusión promedio: ', longitud_difusion_promedio
        WRITE(10, '(A, I0)') '# Número de islas: ', num_islas
        WRITE(10, '(A, I0)') '# Átomos en superficie: ', total_atomos_superficie
        WRITE(10, '(A)') '# Formato: i j altura(i,j)'
        
        DO i = 1, tam_rejilla
            DO j = 1, tam_rejilla
                z_max = find_z_max(material(i,j,:))
                WRITE(10, '(I4, I4, I4)') i, j, z_max
            END DO
        END DO
        CLOSE(10)
    END SUBROUTINE guardar_datos

END PROGRAM simulation_kmc_unificado
