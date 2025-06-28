program kmc_sputtering_simulation
    implicit none

    ! Parámetros
    integer, parameter :: nx = 300  ! Tamaño del sustrato en la dirección x
    integer, parameter :: ny = 300  ! Tamaño del sustrato en la dirección y
    real, parameter :: lambda = 1000.0  ! Tasa promedio de partículas por intervalo (Poisson)
    real, parameter :: kB = 1.380649e-23  ! Constante de Boltzmann en J/K
    real, parameter :: E_diff = 0.6e-19   ! Energía de activación para la difusión de Cr en J
    real, parameter :: E_des = 1.2e-19    ! Energía de activación para la desorción de Cr en J

    ! Variables
    integer :: height(nx, ny)  ! Altura de la película en cada punto del sustrato
    logical :: chemisorbed(nx, ny)  ! Matriz para marcar sitios quimisorbidos
    integer :: i, j, pos_x, pos_y, new_x, new_y
    integer :: nparticles  ! Número de partículas (determinado por Poisson)
    real :: r, energy, temperature, incident_angle
    integer :: seed_size, clock, iostat
    integer, allocatable :: seed(:)

    ! Archivos de salida
    character(len=*), parameter :: output_file = 'height_map_cr.dat'
    character(len=*), parameter :: gnuplot_script = 'plot_script_cr.gnu'

    ! Inicialización
    temperature = 300.0  ! Temperatura del sustrato en Kelvin

    ! Inicializar la semilla aleatoria basada en el reloj del sistema
    ! Esto asegura simulaciones reproducibles pero diferentes en cada ejecución
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    call system_clock(count=clock)
    seed = clock + 37 * [(i - 1, i = 1, seed_size)]
    call random_seed(put=seed)
    deallocate(seed)

    ! Inicializar el sustrato con defectos (5% de probabilidad)
    do j = 1, ny
        do i = 1, nx
            call random_number(r)
            if (r < 0.05) then  ! 5% de probabilidad de defecto
                height(i, j) = 1
                chemisorbed(i, j) = .false.
            else
                height(i, j) = 0
                chemisorbed(i, j) = .false.
            end if
        end do
    end do

    ! Determinar el número de partículas usando la distribución de Poisson
    ! lambda representa la tasa promedio de partículas por intervalo de tiempo
    call random_poisson(lambda, nparticles)
    print *, 'Número de partículas a depositar (Poisson): ', nparticles

    ! Depositar partículas
    do i = 1, nparticles
        ! Generar posiciones aleatorias en el sustrato, asegurando límites
        call random_number(r)
        pos_x = min(max(int(r * nx) + 1, 1), nx)
        call random_number(r)
        pos_y = min(max(int(r * ny) + 1, 1), ny)

        ! Generar ángulo de incidencia aleatorio (0 a π/2 radianes)
        call random_number(r)
        incident_angle = r * 3.14159 / 2.0

        ! Asignar una energía cinética aleatoria a la partícula entrante (0 a 10 eV)
        call random_number(r)
        energy = r * 10.0

        ! Determinar el resultado de la interacción (con ángulo)
        call interact(pos_x, pos_y, energy, height, chemisorbed, nx, ny, incident_angle)

        ! Difusión superficial (con desorción), solo si no está quimisorbido
        if (.not. chemisorbed(pos_x, pos_y)) then
            call diffuse(pos_x, pos_y, height, chemisorbed, nx, ny, temperature)
        end if
    end do

    ! Escribir los resultados en un archivo .dat
    open(unit=10, file=output_file, status='replace', iostat=iostat)
    if (iostat /= 0) then
        print *, 'Error al abrir el archivo ', output_file
        stop
    end if
    do j = 1, ny
        do i = 1, nx
            write(10, *) i, j, height(i, j)
        end do
        write(10, *)  ! Nueva línea para separar filas
    end do
    close(10)

    ! Generar script para gnuplot
    open(unit=20, file=gnuplot_script, status='replace', iostat=iostat)
    if (iostat /= 0) then
        print *, 'Error al abrir el archivo ', gnuplot_script
        stop
    end if
    write(20, *) 'set terminal png'
    write(20, *) 'set output "height_map_cr.png"'
    write(20, *) 'set title "Height Map after ', nparticles, ' particles (Cr)"'
    write(20, *) 'set xlabel "X"'
    write(20, *) 'set ylabel "Y"'
    write(20, *) 'plot "', trim(output_file), '" using 1:2:3 with image'
    close(20)

    print *, 'Simulación completada. Los resultados se han guardado en ', output_file
    print *, 'Script de gnuplot creado. Ejecuta gnuplot para visualizar los resultados.'

contains

    ! Subrutina para generar un número aleatorio de Poisson (algoritmo de Knuth)
    ! Genera el número de partículas basado en la tasa promedio lambda
    subroutine random_poisson(lambda, k)
        real, intent(in) :: lambda
        integer, intent(out) :: k
        real :: L, p, u
        L = exp(-lambda)
        k = 0
        p = 1.0
        do
            call random_number(u)
            p = p * u
            if (p < L) exit
            k = k + 1
        end do
    end subroutine random_poisson

    ! Subrutina para interactuar con el sustrato
    subroutine interact(x, y, energy, height, chemisorbed, nx, ny, angle)
        integer, intent(inout) :: x, y
        real, intent(in) :: energy, angle
        integer, intent(inout) :: height(nx, ny)
        logical, intent(inout) :: chemisorbed(nx, ny)
        integer, intent(in) :: nx, ny
        real :: adsorption_prob, reflection_prob, chemisorption_prob, r, total_prob
        real, parameter :: E_ads = 1.2e-19  ! Energía de activación para la adsorción de Cr en J
        real, parameter :: E_sput = 2.2e-19 ! Energía de activación para el sputtering de Cr en J

        ! Calcular probabilidades basadas en energía, ángulo y propiedades del Cr
        ! La probabilidad de adsorción incluye un factor de Boltzmann para E_ads
        adsorption_prob = 0.7 * exp(-E_ads / (kB * 300.0)) * (1.0 - sin(angle) * (energy / 10.0))
        reflection_prob = 0.2
        chemisorption_prob = 0.1 * exp(-E_ads / (kB * 300.0))  ! Ajustado para Cr
        total_prob = adsorption_prob + reflection_prob + chemisorption_prob

        ! Normalizar probabilidades para que sumen 1
        adsorption_prob = adsorption_prob / total_prob
        reflection_prob = reflection_prob / total_prob
        chemisorption_prob = chemisorption_prob / total_prob

        call random_number(r)
        if (r < adsorption_prob) then
            ! Adsorción: incrementar altura
            height(x, y) = height(x, y) + 1
            chemisorbed(x, y) = .false.
        else if (r < adsorption_prob + reflection_prob) then
            ! Reflexión: no hacer nada
        else if (r < adsorption_prob + reflection_prob + chemisorption_prob) then
            ! Quimisorción: incrementar altura y marcar como quimisorbido
            height(x, y) = height(x, y) + 1
            chemisorbed(x, y) = .true.
        else
            ! Sputtering: remover material si no está quimisorbido
            if (.not. chemisorbed(x, y)) then
                call sputter(x, y, height, chemisorbed, nx, ny, angle)
            end if
        end if
    end subroutine interact

    ! Subrutina para sputtering
    subroutine sputter(x, y, height, chemisorbed, nx, ny, angle)
        integer, intent(in) :: x, y
        integer, intent(inout) :: height(nx, ny)
        logical, intent(in) :: chemisorbed(nx, ny)
        integer, intent(in) :: nx, ny
        real, intent(in) :: angle
        integer :: dx, dy, attempts = 0
        real :: r

        ! Intentar hasta 5 veces encontrar un sitio para sputtering
        do while (attempts < 5)
            attempts = attempts + 1
            call random_number(r)
            dx = int(cos(angle) * 3) - 1  ! Dirección basada en el ángulo
            call random_number(r)
            dy = int(sin(angle) * 3) - 1

            ! Aplicar condiciones de contorno periódicas
            new_x = modulo(x + dx - 1, nx) + 1
            new_y = modulo(y + dy - 1, ny) + 1

            ! Verificar si el sitio tiene material y no está quimisorbido
            if (height(new_x, new_y) > 0 .and. .not. chemisorbed(new_x, new_y)) then
                height(new_x, new_y) = height(new_x, new_y) - 1
                exit
            end if
        end do
    end subroutine sputter

    ! Subrutina para difusión superficial
    subroutine diffuse(x, y, height, chemisorbed, nx, ny, temperature)
        integer, intent(inout) :: x, y
        integer, intent(inout) :: height(nx, ny)
        logical, intent(in) :: chemisorbed(nx, ny)
        integer, intent(in) :: nx, ny
        real, intent(in) :: temperature
        integer :: dx, dy
        real :: r, diffusion_probability, energy_diff, desorption_prob
        real :: boltzmann_factor

        ! Solo permitir difusión si el sitio no está quimisorbido
        if (chemisorbed(x, y)) return

        ! Intentar difundir a un sitio vecino
        call random_number(r)
        dx = int(r * 3) - 1  ! -1, 0, 1
        call random_number(r)
        dy = int(r * 3) - 1

        ! Aplicar condiciones de contorno periódicas
        new_x = modulo(x + dx - 1, nx) + 1
        new_y = modulo(y + dy - 1, ny) + 1

        ! Calcular probabilidad de difusión usando el factor de Boltzmann
        energy_diff = E_diff * (height(x, y) - height(new_x, new_y))
        boltzmann_factor = exp(-energy_diff / (kB * temperature))
        diffusion_probability = min(1.0, boltzmann_factor)

        call random_number(r)
        if (height(new_x, new_y) < height(x, y) .and. r < diffusion_probability .and. .not. chemisorbed(new_x, new_y)) then
            ! Difusión: mover partícula al sitio vecino
            height(x, y) = height(x, y) - 1
            height(new_x, new_y) = height(new_x, new_y) + 1
            x = new_x
            y = new_y
        else
            ! Calcular probabilidad de desorción usando E_des
            desorption_prob = exp(-E_des / (kB * temperature))
            call random_number(r)
            if (r < desorption_prob) then
                ! Desorción: remover partícula
                height(x, y) = height(x, y) - 1
            end if
        end if
    end subroutine diffuse

end program kmc_sputtering_simulation
