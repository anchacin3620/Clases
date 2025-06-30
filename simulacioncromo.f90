program MonteCarloSputtering3D_Cr
    implicit none

    ! Parámetros ajustables del sistema
    integer, parameter :: nx = 300        ! Tamaño del sustrato en x
    integer, parameter :: ny = 300        ! Tamaño del sustrato en y
    integer, parameter :: Lz = 100        ! Altura máxima en z
    integer, parameter :: Mg_thickness = 2 ! Espesor del sustrato Mg (capas)
    real, parameter :: lambda = 20000.0   ! Tasa promedio de partículas por intervalo (Poisson)
    real, parameter :: PI = 3.1415926535  ! Constante pi
    real, parameter :: kB = 1.380649e-23  ! Constante de Boltzmann (J/K)
    real, parameter :: eV_to_J = 1.602e-19 ! Conversión eV a J

    ! Parámetros físicos ajustables (específicos para Cr sobre Mg)
    real, parameter :: T_subs = 500.0     ! Temperatura del sustrato (K)
    real, parameter :: T_plasma = 1000.0  ! Temperatura del plasma en sputtering (K)
    real, parameter :: E_d_Cr = 0.6       ! Barrera de energía para difusión Cr/Mg (eV)
    real, parameter :: E_des_Cr = 1.2     ! Energía de desorción Cr/Mg (eV)
    real, parameter :: E_ads_Cr = 1.2     ! Energía de adsorción Cr/Mg (eV)
    real, parameter :: E_sput_Cr = 2.2    ! Energía de sputtering Cr/Mg (eV)
    real, parameter :: sticking_Cr = 0.80 ! Coeficiente de adherencia Cr sobre Mg

    ! Variables de la simulación
    integer :: substrate(nx, ny, Lz)      ! Matriz 3D (0 = vacío, 2 = Cr, 3 = Mg)
    logical :: chemisorbed(nx, ny, Lz)    ! Matriz para sitios quimisorbidos
    integer :: height(nx, ny)             ! Altura de la película en cada (x, y)
    real :: theta, phi                    ! Ángulos de incidencia (radianes)
    real :: x_pos, y_pos, energy          ! Posición y energía de llegada
    integer :: i, j, k, l, xpos, ypos, zpos ! Índices y posiciones
    real :: r1, r2, r3, r4                ! Números aleatorios
    real :: diff_prob_Cr                  ! Probabilidad de difusión para Cr
    integer :: nparticles                 ! Número de partículas (Poisson)
    integer :: seed_size, clock, iostat   ! Para semilla aleatoria y manejo de errores
    integer, allocatable :: seed(:)

    ! Inicialización
    ! Inicializar semilla aleatoria con reloj del sistema para reproducibilidad
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    call system_clock(count=clock)
    seed = clock + 37 * [(i - 1, i = 1, seed_size)]
    call random_seed(put=seed)
    deallocate(seed)

    ! Inicializar matrices
    substrate = 0                         ! Espacio vacío inicialmente
    chemisorbed = .false.                 ! Sin sitios quimisorbidos inicialmente

    ! Inicializar sustrato Mg (001) con espesor Mg_thickness
    do i = 1, nx
        do j = 1, ny
            do k = 1, Mg_thickness
                substrate(i, j, k) = 3    ! Mg = 3
                chemisorbed(i, j, k) = .false.
            end do
            height(i, j) = Mg_thickness   ! Altura inicial
        end do
    end do

    ! Probabilidad de difusión para Cr (Arrhenius)
    diff_prob_Cr = exp(-E_d_Cr * eV_to_J / (kB * T_subs))

    ! Determinar número de partículas con distribución de Poisson
    call random_poisson(lambda, nparticles)
    print *, 'Número de partículas a depositar (Poisson): ', nparticles

    ! Archivos de salida
    open(unit=10, file='Cr_height.dat', status='replace', iostat=iostat)
    if (iostat /= 0) then
        print *, 'Error al abrir el archivo Cr_height.dat'
        stop
    end if
    open(unit=11, file='Cr_3d.vtk', status='replace', iostat=iostat)
    if (iostat /= 0) then
        print *, 'Error al abrir el archivo Cr_3d.vtk'
        stop
    end if

    ! Simulación Monte Carlo Cinético 3D para Cr
    do k = 1, nparticles
        ! 1. Generar posición y ángulos (sputtering)
        call random_number(r1)
        call random_number(r2)
        x_pos = r1 * real(nx)
        y_pos = r2 * real(ny)
        xpos = min(max(int(x_pos) + 1, 1), nx)
        ypos = min(max(int(y_pos) + 1, 1), ny)

        call random_number(r3)
        theta = acos(1.0 - 2.0 * r3)       ! Distribución coseno para theta
        call random_number(r4)
        phi = 2.0 * PI * r4                ! Ángulo azimutal uniforme

        ! 2. Energía según Maxwell-Boltzmann
        call random_number(r1)
        energy = -kB * T_plasma / eV_to_J * log(r1)  ! Distribución Boltzmann (eV)

        ! 3. Determinar interacción (adsorción, reflexión, quimisorción, sputtering)
        call interact(xpos, ypos, energy, theta, substrate, chemisorbed, height, nx, ny, Lz)

        ! 4. Progreso
        if (mod(k, 1000) == 0) then
            print *, 'Cr depositadas:', k
        end if
    end do

    ! Escribir perfil de altura proyectado (Cr_height.dat)
    write(10, *) '# x  y  height  material'
    do i = 1, nx
        do j = 1, ny
            write(10, '(I5, I5, I5, I5)') i, j, height(i, j), substrate(i, j, height(i, j))
        end do
        write(10, *)  ! Línea en blanco para separar bloques en Gnuplot
    end do

    ! Escribir datos 3D en formato VTK simplificado (Cr_3d.vtk)
    write(11, '(A)') '# vtk DataFile Version 3.0'
    write(11, '(A)') 'Sputtering 3D Cr on Mg(001) 300x300'
    write(11, '(A)') 'ASCII'
    write(11, '(A)') 'DATASET STRUCTURED_POINTS'
    write(11, '(A, I4, I4, I4)') 'DIMENSIONS ', nx, ny, Lz
    write(11, '(A)') 'ORIGIN 0 0 0'
    write(11, '(A)') 'SPACING 1 1 1'
    write(11, '(A, I10)') 'POINT_DATA ', nx * ny * Lz
    write(11, '(A)') 'SCALARS material int 1'
    write(11, '(A)') 'LOOKUP_TABLE default'
    do l = 1, Lz
        do j = 1, ny
            do i = 1, nx
                write(11, '(I2)') substrate(i, j, l)
            end do
        end do
    end do

    ! Estadísticas
    print *, 'Altura promedio Cr sobre Mg:', real(sum(height) - Mg_thickness * nx * ny) / real(nx * ny)
    print *, 'Probabilidad de difusión Cr:', diff_prob_Cr
    print *, 'Simulación completada. Datos en Cr_height.dat y Cr_3d.vtk'

    ! Cerrar archivos
    close(10)
    close(11)

contains

    ! Subrutina para generar un número aleatorio de Poisson (algoritmo de Knuth)
    ! Genera el número de partículas basado en la tasa promedio lambda
    subroutine random_poisson(lambda, k)
        real, intent(in) :: lambda
        integer, intent(out) :: k
        real :: L, p, u
        integer :: iter
        parameter (iter = 20000)  ! Límite máximo de iteraciones
        L = exp(-lambda)
        k = 0
        p = 1.0
        do i = 1, iter
            call random_number(u)
            p = p * u
            if (p < L) exit
            k = k + 1
            if (i == iter) then
                print *, 'Warning: random_poisson reached iteration limit, k set to ', k
                exit
            end if
        end do
    end subroutine random_poisson

    ! Subrutina para determinar la interacción de la partícula con el sustrato
    subroutine interact(xpos, ypos, energy, theta, substrate, chemisorbed, height, nx, ny, Lz)
        integer, intent(in) :: nx, ny, Lz
        integer, intent(inout) :: xpos, ypos
        real, intent(in) :: energy, theta
        integer, intent(inout) :: substrate(nx, ny, Lz)
        logical, intent(inout) :: chemisorbed(nx, ny, Lz)
        integer, intent(inout) :: height(nx, ny)
        real :: adsorption_prob, reflection_prob, chemisorption_prob, sputter_prob, total_prob, r
        integer :: i, j, zpos

        ! Calcular probabilidades basadas en energía, ángulo y propiedades del Cr
        adsorption_prob = sticking_Cr * exp(-E_ads_Cr * eV_to_J / (kB * T_plasma)) * (1.0 - sin(theta))
        reflection_prob = 0.2 * (1.0 - exp(-energy / (kB * T_plasma / eV_to_J)))
        chemisorption_prob = 0.1 * exp(-E_ads_Cr * eV_to_J / (kB * T_plasma))
        sputter_prob = 0.05 * (energy / (E_sput_Cr * eV_to_J))  ! Depende de la energía
        total_prob = adsorption_prob + reflection_prob + chemisorption_prob + sputter_prob

        ! Normalizar probabilidades
        adsorption_prob = adsorption_prob / total_prob
        reflection_prob = reflection_prob / total_prob
        chemisorption_prob = chemisorption_prob / total_prob
        sputter_prob = sputter_prob / total_prob

        call random_number(r)
        if (r < adsorption_prob) then
            ! Adsorción: depositar Cr en el punto más alto
            zpos = height(xpos, ypos) + 1
            do i = max(1, xpos-1), min(nx, xpos+1)
                do j = max(1, ypos-1), min(ny, ypos+1)
                    if (height(i, j) >= zpos) then
                        zpos = height(i, j) + 1  ! Ajustar al pico más alto
                    end if
                end do
            end do
            if (zpos <= Lz) then
                substrate(xpos, ypos, zpos) = 2  ! Cr = 2
                chemisorbed(xpos, ypos, zpos) = .false.
                height(xpos, ypos) = max(height(xpos, ypos), zpos)
                ! Difusión superficial si está por encima del sustrato Mg
                if (zpos > Mg_thickness) then
                    call diffuse(xpos, ypos, zpos, substrate, chemisorbed, height, nx, ny, Lz)
                end if
            end if
        else if (r < adsorption_prob + reflection_prob) then
            ! Reflexión: no hacer nada
        else if (r < adsorption_prob + reflection_prob + chemisorption_prob) then
            ! Quimisorción: depositar Cr y marcar como quimisorbido
            zpos = height(xpos, ypos) + 1
            do i = max(1, xpos-1), min(nx, xpos+1)
                do j = max(1, ypos-1), min(ny, ypos+1)
                    if (height(i, j) >= zpos) then
                        zpos = height(i, j) + 1
                    end if
                end do
            end do
            if (zpos <= Lz) then
                substrate(xpos, ypos, zpos) = 2  ! Cr = 2
                chemisorbed(xpos, ypos, zpos) = .true.
                height(xpos, ypos) = max(height(xpos, ypos), zpos)
            end if
        else
            ! Sputtering: remover material si no está quimisorbido
            call sputter(xpos, ypos, substrate, chemisorbed, height, nx, ny, Lz, theta)
        end if
    end subroutine interact

    ! Subrutina para difusión superficial en 3D
    subroutine diffuse(xpos, ypos, zpos, substrate, chemisorbed, height, nx, ny, Lz)
        integer, intent(in) :: nx, ny, Lz
        integer, intent(inout) :: xpos, ypos, zpos
        integer, intent(inout) :: substrate(nx, ny, Lz)
        logical, intent(in) :: chemisorbed(nx, ny, Lz)
        integer, intent(inout) :: height(nx, ny)
        real :: r1, r2, desorption_prob
        integer :: new_x, new_y

        ! Solo permitir difusión si el sitio no está quimisorbido
        if (chemisorbed(xpos, ypos, zpos)) return

        call random_number(r1)
        if (r1 < diff_prob_Cr) then
            call random_number(r2)
            if (r2 < 0.25 .and. height(modulo(xpos-2, nx)+1, ypos) < zpos) then
                new_x = modulo(xpos-2, nx) + 1  ! Mover a la izquierda (periódico)
                substrate(xpos, ypos, zpos) = 0
                substrate(new_x, ypos, zpos) = 2
                xpos = new_x
            else if (r2 < 0.5 .and. height(modulo(xpos, nx)+1, ypos) < zpos) then
                new_x = modulo(xpos, nx) + 1    ! Mover a la derecha
                substrate(xpos, ypos, zpos) = 0
                substrate(new_x, ypos, zpos) = 2
                xpos = new_x
            else if (r2 < 0.75 .and. height(xpos, modulo(ypos-2, ny)+1) < zpos) then
                new_y = modulo(ypos-2, ny) + 1  ! Mover hacia atrás
                substrate(xpos, ypos, zpos) = 0
                substrate(xpos, new_y, zpos) = 2
                ypos = new_y
            else if (height(xpos, modulo(ypos, ny)+1) < zpos) then
                new_y = modulo(ypos, ny) + 1    ! Mover hacia adelante
                substrate(xpos, ypos, zpos) = 0
                substrate(xpos, new_y, zpos) = 2
                ypos = new_y
            end if
            height(xpos, ypos) = max(height(xpos, ypos), zpos)  ! Actualizar altura
        else
            ! Calcular probabilidad de desorción
            desorption_prob = exp(-E_des_Cr * eV_to_J / (kB * T_subs))
            call random_number(r1)
            if (r1 < desorption_prob) then
                substrate(xpos, ypos, zpos) = 0  ! Desorción: remover partícula
                if (zpos == height(xpos, ypos)) then
                    ! Recalcular altura si se remueve la partícula superior
                    do k = zpos-1, 1, -1
                        if (substrate(xpos, ypos, k) /= 0) then
                            height(xpos, ypos) = k
                            exit
                        end if
                    end do
                    if (k == 0) height(xpos, ypos) = Mg_thickness
                end if
            end if
        end if
    end subroutine diffuse

    ! Subrutina para sputtering
    subroutine sputter(xpos, ypos, substrate, chemisorbed, height, nx, ny, Lz, theta)
        integer, intent(in) :: nx, ny, Lz
        integer, intent(in) :: xpos, ypos
        integer, intent(inout) :: substrate(nx, ny, Lz)
        logical, intent(in) :: chemisorbed(nx, ny, Lz)
        integer, intent(inout) :: height(nx, ny)
        real, intent(in) :: theta
        integer :: new_x, new_y, zpos, attempts
        real :: r

        ! Intentar hasta 5 veces encontrar un sitio para sputtering
        attempts = 0
        zpos = height(xpos, ypos)
        do while (attempts < 5)
            attempts = attempts + 1
            call random_number(r)
            new_x = modulo(xpos + int(cos(theta) * 3) - 1, nx) + 1
            call random_number(r)
            new_y = modulo(ypos + int(sin(theta) * 3) - 1, ny) + 1
            if (height(new_x, new_y) > Mg_thickness .and. .not. chemisorbed(new_x, new_y, height(new_x, new_y))) then
                substrate(new_x, new_y, height(new_x, new_y)) = 0
                ! Recalcular altura
                do k = height(new_x, new_y)-1, 1, -1
                    if (substrate(new_x, new_y, k) /= 0) then
                        height(new_x, new_y) = k
                        exit
                    end if
                end do
                if (k == 0) height(new_x, new_y) = Mg_thickness
                exit
            end if
        end do
    end subroutine sputter

end program MonteCarloSputtering3D_Cr
