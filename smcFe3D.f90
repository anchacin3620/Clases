program MonteCarloSputtering3D_Fe
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

    ! Parámetros físicos ajustables (específicos para Fe sobre Mg)
    real, parameter :: T_subs = 500.0     ! Temperatura del sustrato (K)
    real, parameter :: T_plasma = 1000.0  ! Temperatura del plasma en sputtering (K)
    real, parameter :: E_d_Fe = 0.7       ! Barrera de energía para difusión Fe/Mg (eV)
    real, parameter :: E_des_Fe = 1.4     ! Energía de desorción Fe/Mg (eV)
    real, parameter :: E_ads_Fe = 1.3     ! Energía de adsorción Fe/Mg (eV)
    real, parameter :: E_sput_Fe = 2.3    ! Energía de sputtering Fe/Mg (eV)
    real, parameter :: sticking_Fe = 0.85 ! Coeficiente de adherencia Fe sobre Mg

    ! Variables de la simulación
    integer :: substrate(nx, ny, Lz)      ! Matriz 3D (0 = vacío, 2 = Fe, 3 = Mg)
    logical :: chemisorbed(nx, ny, Lz)    ! Matriz para sitios quimisorbidos
    integer :: height(nx, ny)             ! Altura de la película en cada (x, y)
    real :: theta, phi                    ! Ángulos de incidencia (radianes)
    real :: x_pos, y_pos, energy          ! Posición y energía de llegada
    integer :: i, j, k, l, xpos, ypos, zpos ! Índices y posiciones
    real :: r1, r2, r3, r4                ! Números aleatorios
    real :: diff_prob_Fe                  ! Probabilidad de difusión para Fe
    integer :: nparticles                 ! Número de partículas (Poisson)
    integer :: seed_size, clock, iostat   ! Para semilla aleatoria y manejo de errores
    integer, allocatable :: seed(:)

    ! Inicialización
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    call system_clock(count=clock)
    seed = clock + 37 * [(i - 1, i = 1, seed_size)]
    call random_seed(put=seed)
    deallocate(seed)

    substrate = 0
    chemisorbed = .false.

    do i = 1, nx
        do j = 1, ny
            do k = 1, Mg_thickness
                substrate(i, j, k) = 3
                chemisorbed(i, j, k) = .false.
            end do
            height(i, j) = Mg_thickness
        end do
    end do

    diff_prob_Fe = exp(-E_d_Fe * eV_to_J / (kB * T_subs))

    call random_poisson(lambda, nparticles)
    print *, 'Número de partículas a depositar (Poisson): ', nparticles

    open(unit=10, file='Fe_height.dat', status='replace', iostat=iostat)
    if (iostat /= 0) then
        print *, 'Error al abrir el archivo Fe_height.dat'
        stop
    end if
    open(unit=11, file='Fe_3d_slice_z10.dat', status='replace', iostat=iostat) ! Ejemplo de slice en z=10
    if (iostat /= 0) then
        print *, 'Error al abrir el archivo Fe_3d_slice_z10.dat'
        stop
    end if

    do k = 1, nparticles
        call random_number(r1)
        call random_number(r2)
        x_pos = r1 * real(nx)
        y_pos = r2 * real(ny)
        xpos = min(max(int(x_pos) + 1, 1), nx)
        ypos = min(max(int(y_pos) + 1, 1), ny)

        call random_number(r3)
        theta = acos(1.0 - 2.0 * r3)
        call random_number(r4)
        phi = 2.0 * PI * r4

        call random_number(r1)
        energy = -kB * T_plasma / eV_to_J * log(r1)

        call interact(xpos, ypos, energy, theta, substrate, chemisorbed, height, nx, ny, Lz)

        if (mod(k, 1000) == 0) then
            print *, 'Fe depositadas:', k
        end if
    end do

    write(10, *) '# x  y  height  material'
    do i = 1, nx
        do j = 1, ny
            write(10, '(I5, I5, I5, I5)') i, j, height(i, j), substrate(i, j, height(i, j))
        end do
        write(10, *)
    end do

    ! Escribir una slice 2D en z=10 como ejemplo
    write(11, *) '# x  y  material'
    do i = 1, nx
        do j = 1, ny
            write(11, '(I5, I5, I5)') i, j, substrate(i, j, 10)
        end do
        write(11, *)
    end do

    print *, 'Altura promedio Fe sobre Mg:', real(sum(height) - Mg_thickness * nx * ny) / real(nx * ny)
    print *, 'Probabilidad de difusión Fe:', diff_prob_Fe
    print *, 'Simulación completada. Datos en Fe_height.dat y Fe_3d_slice_z10.dat'

    close(10)
    close(11)

contains

    subroutine random_poisson(lambda, k)
        real, intent(in) :: lambda
        integer, intent(out) :: k
        real :: L, p, u
        integer :: iter
        parameter (iter = 20000)
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

    subroutine interact(xpos, ypos, energy, theta, substrate, chemisorbed, height, nx, ny, Lz)
        integer, intent(in) :: nx, ny, Lz
        integer, intent(inout) :: xpos, ypos
        real, intent(in) :: energy, theta
        integer, intent(inout) :: substrate(nx, ny, Lz)
        logical, intent(inout) :: chemisorbed(nx, ny, Lz)
        integer, intent(inout) :: height(nx, ny)
        real :: adsorption_prob, reflection_prob, chemisorption_prob, sputter_prob, total_prob, r
        integer :: i, j, zpos

        adsorption_prob = sticking_Fe * exp(-E_ads_Fe * eV_to_J / (kB * T_plasma)) * (1.0 - sin(theta))
        reflection_prob = 0.2 * (1.0 - exp(-energy / (kB * T_plasma / eV_to_J)))
        chemisorption_prob = 0.1 * exp(-E_ads_Fe * eV_to_J / (kB * T_plasma))
        sputter_prob = 0.05 * (energy / (E_sput_Fe * eV_to_J))
        total_prob = adsorption_prob + reflection_prob + chemisorption_prob + sputter_prob

        adsorption_prob = adsorption_prob / total_prob
        reflection_prob = reflection_prob / total_prob
        chemisorption_prob = chemisorption_prob / total_prob
        sputter_prob = sputter_prob / total_prob

        call random_number(r)
        if (r < adsorption_prob) then
            zpos = height(xpos, ypos) + 1
            do i = max(1, xpos-1), min(nx, xpos+1)
                do j = max(1, ypos-1), min(ny, ypos+1)
                    if (height(i, j) >= zpos) zpos = height(i, j) + 1
                end do
            end do
            if (zpos <= Lz) then
                substrate(xpos, ypos, zpos) = 2
                chemisorbed(xpos, ypos, zpos) = .false.
                height(xpos, ypos) = max(height(xpos, ypos), zpos)
                if (zpos > Mg_thickness) call diffuse(xpos, ypos, zpos, substrate, chemisorbed, height, nx, ny, Lz)
            end if
        else if (r < adsorption_prob + reflection_prob) then
        else if (r < adsorption_prob + reflection_prob + chemisorption_prob) then
            zpos = height(xpos, ypos) + 1
            do i = max(1, xpos-1), min(nx, xpos+1)
                do j = max(1, ypos-1), min(ny, ypos+1)
                    if (height(i, j) >= zpos) zpos = height(i, j) + 1
                end do
            end do
            if (zpos <= Lz) then
                substrate(xpos, ypos, zpos) = 2
                chemisorbed(xpos, ypos, zpos) = .true.
                height(xpos, ypos) = max(height(xpos, ypos), zpos)
            end if
        else
            call sputter(xpos, ypos, substrate, chemisorbed, height, nx, ny, Lz, theta)
        end if
    end subroutine interact

    subroutine diffuse(xpos, ypos, zpos, substrate, chemisorbed, height, nx, ny, Lz)
        integer, intent(in) :: nx, ny, Lz
        integer, intent(inout) :: xpos, ypos, zpos
        integer, intent(inout) :: substrate(nx, ny, Lz)
        logical, intent(in) :: chemisorbed(nx, ny, Lz)
        integer, intent(inout) :: height(nx, ny)
        real :: r1, r2, desorption_prob
        integer :: new_x, new_y

        if (chemisorbed(xpos, ypos, zpos)) return

        call random_number(r1)
        if (r1 < diff_prob_Fe) then
            call random_number(r2)
            if (r2 < 0.25 .and. height(modulo(xpos-2, nx)+1, ypos) < zpos) then
                new_x = modulo(xpos-2, nx) + 1
                substrate(xpos, ypos, zpos) = 0
                substrate(new_x, ypos, zpos) = 2
                xpos = new_x
            else if (r2 < 0.5 .and. height(modulo(xpos, nx)+1, ypos) < zpos) then
                new_x = modulo(xpos, nx) + 1
                substrate(xpos, ypos, zpos) = 0
                substrate(new_x, ypos, zpos) = 2
                xpos = new_x
            else if (r2 < 0.75 .and. height(xpos, modulo(ypos-2, ny)+1) < zpos) then
                new_y = modulo(ypos-2, ny) + 1
                substrate(xpos, ypos, zpos) = 0
                substrate(xpos, new_y, zpos) = 2
                ypos = new_y
            else if (height(xpos, modulo(ypos, ny)+1) < zpos) then
                new_y = modulo(ypos, ny) + 1
                substrate(xpos, ypos, zpos) = 0
                substrate(xpos, new_y, zpos) = 2
                ypos = new_y
            end if
            height(xpos, ypos) = max(height(xpos, ypos), zpos)
        else
            desorption_prob = exp(-E_des_Fe * eV_to_J / (kB * T_subs))
            call random_number(r1)
            if (r1 < desorption_prob) then
                substrate(xpos, ypos, zpos) = 0
                if (zpos == height(xpos, ypos)) then
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

    subroutine sputter(xpos, ypos, substrate, chemisorbed, height, nx, ny, Lz, theta)
        integer, intent(in) :: nx, ny, Lz
        integer, intent(in) :: xpos, ypos
        integer, intent(inout) :: substrate(nx, ny, Lz)
        logical, intent(in) :: chemisorbed(nx, ny, Lz)
        integer, intent(inout) :: height(nx, ny)
        real, intent(in) :: theta
        integer :: new_x, new_y, zpos, attempts
        real :: r

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

end program MonteCarloSputtering3D_Fe

! Gnuplot scripts for visualization
! 2D Height Map Plot
! Save as 'plot_2d.gp' and run with 'gnuplot plot_2d.gp'
! set datafile separator whitespace
! set title "2D Height Map of Fe on Mg Substrate"
! set xlabel "X"
! set ylabel "Y"
! set zlabel "Height"
! set grid
! plot 'Fe_height.dat' using 1:2:3 with linespoints
! pause -1 "Press Enter to exit"

! 3D Surface Plot
! Save as 'plot_3d.gp' and run with 'gnuplot plot_3d.gp'
! set datafile separator whitespace
! set title "3D Surface of Fe on Mg Substrate"
! set xlabel "X"
! set ylabel "Y"
! set zlabel "Height"
! set grid
! set pm3d map
! splot 'Fe_height.dat' using 1:2:3 with pm3d
! pause -1 "Press Enter to exit"

! 3D Slice Plot (e.g., z=10)
! Save as 'plot_slice.gp' and run with 'gnuplot plot_slice.gp'
! set datafile separator whitespace
! set title "Slice at Z=10 of Fe on Mg Substrate"
! set xlabel "X"
! set ylabel "Y"
! set zlabel "Material"
! set pm3d map
! splot 'Fe_3d_slice_z10.dat' using 1:2:3 with pm3d
! pause -1 "Press Enter to exit"
