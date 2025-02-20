! Serial version
program jacobi_serial
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    integer :: i, j, k, istat
    integer :: Nx, Ny, Lx, Ly, dx, dy, x1, x2, grid_area
    integer :: nIterations, nIterationsOutput 
    
    real(kind=dp) :: A, alpha, r2
    real(kind=dp) :: potential, free_parameter, tolerance
    real(kind=dp), allocatable, dimension(:,:) :: grid, previous_grid, charge_density
    
    logical :: all_converged
    real(kind=dp) :: start_time, finish_time 

    ! Define grid: 
    Nx = 500                        ! Number of grid points in x
    Ny = 500                        ! Number of grid points in y
    Lx = 1.0                        ! Length of domain in x 
    Ly = 1.0                        ! Length of domain in x 
    dx = Lx/(Nx-1)                  ! grid spacing in x 
    dy = Ly/(Ny-1)                  ! grid spacing in y 
    grid_area = Nx*Ny        
    x1 = (Nx+1)/2 
    x2 = (Ny+1)/2

    allocate(grid(0:Nx+1,0:Ny+1),stat=istat)
    if (istat.ne.0) stop "Error allocating grid array"
    allocate(previous_grid(0:Nx+1,0:Ny+1),stat=istat)
    if (istat.ne.0) stop "Error allocating previous_grid array"
    allocate(charge_density(0:Nx+1,0:Ny+1),stat=istat)
    if (istat.ne.0) stop "Error allocating charge_density array"

    grid(:,:) = 0.0_dp 
    previous_grid(:,:) = 0.0_dp 
    charge_density(:,:) = 0.0_dp 

    ! Simulation parameters: 
    nIterations = 100000                ! Number of time steps
    nIterationsOutput = nIterations/5   ! Number of times we output a file
    free_parameter = 0.8_dp             ! Overelaxation factor (must be less than one for overrelaxation)
    tolerance = 1.0E-8_dp               ! Tolerance to accept convergence
    start_time = 0.0_dp 
    finish_time = 0.0_dp
    
    ! Gaussian distribution 
    A = 0.1_dp                      ! Max value of potential
    alpha = 0.01_dp                 ! Width of gaussian 

    do i = 1, Nx
        do j = 1, Ny
            r2 = (real(i-x1, dp)**2 + real(j-x2, dp)**2)
            ! Initialise the field with either 0 or gaussian -> This helps for a smooth convergence.
            grid(i,j) = A * exp(-alpha * r2)
            ! Initialise the charge density with a uniform or gaussian distribution. 
            charge_density(i,j) = A * (2.0_dp * alpha - 4.0_dp * alpha**2 * r2) * exp(-alpha * r2)
            ! charge_density(i,j) = 3.0_dp 
        end do
    end do

    ! Print an image of the initial conditions
    call write_to_output(0)
    call cpu_time(start_time)

    convergence_loop: do k = 1, nIterations 
        ! Store the previous grid
        previous_grid(:,:) = grid(:,:)

        ! Jacobi: 

        ! The grid is initialised with 0 values. By iterating between 1 and Nx/Ny - 1, we are maintaining a boundary of 0, 
        ! therefore maintaining Dirichlet boundary conditions. 
        do j = 2, Nx-1
            do i = 2, Ny-1  
                ! Now update the potential value using the nearest neighbours to grid(i,j) and the charge density
                potential = 0.25_dp*(grid(i+1,j) + grid(i-1,j) + grid(i,j+1) + grid(i,j-1) + (dx**2) * charge_density(i,j))
                ! Update grid
                grid(i,j) = grid(i,j) + free_parameter*(potential - grid(i,j))
            end do 
        end do     
    
        ! Periodically write to output 
        if (mod(k,nIterationsOutput).eq.0) then 
            call write_to_output(k)
        end if 
        
        
        call check_convergence 
        
        if (all_converged) print*, "Converged in", k, "iterations"
        if (all_converged) call write_to_output(k)
        if (all_converged) exit convergence_loop
    
    end do convergence_loop

    ! Print an image of the final conditions
    call cpu_time(finish_time)

    print*, "Time to compute:", finish_time-start_time, "s"


    deallocate(grid,stat=istat)
    if (istat.ne.0) stop 'Error deallocating grid array'
    deallocate(previous_grid,stat=istat)
    if (istat.ne.0) stop 'Error deallocating previous_grid array'
    deallocate(charge_density,stat=istat)
    if (istat.ne.0) stop 'Error deallocating previous_grid array'
  
contains 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to write to output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_to_output(k)
    implicit none 
    integer :: k
    character(len=20) :: filename 
    character(len=24) :: filename_pgm, filename_dat 
    character(len=5) :: k_str

    ! Convert integer k to string
    write(k_str, '(I0)') k  ! 'I0' format writes integer without leading spaces

    ! Construct the output filenames based on k
    filename = "phi_" // trim(k_str)
    ! out_unit = k 
    filename_pgm = trim(filename) // ".pgm"
    filename_dat = trim(filename) // ".dat"
            
    ! Visualise 2D domain 
    call create_pgm(Nx,Ny,grid,filename_pgm)

    open(file=filename_dat,unit=k,status='unknown',iostat=istat)
    if (istat.ne.0) stop "Error opening write_to_output .dat"
    write(k,fmt=*,iostat=istat) grid(:,:)
    if (istat.ne.0) stop "Error writing to write_to_output .dat"
    close(unit=k,iostat=istat)
    if (istat.ne.0) stop "Error closing write_to_output .dat"
    end subroutine write_to_output 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the grid has converged
    ! by comparing the previous value point at the current grid point
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine check_convergence
            implicit none
            integer,parameter :: dp = selected_real_kind(15,300)
            integer :: converged 
            
            converged = 0
            all_converged = .false.
            ! This is an expensive test
            do i = 1, Nx
                do j = 1, Ny
                    ! Check that each point value is within some tolerance of the previous value of that point
                    if (abs(grid(i,j)-previous_grid(i,j)) .lt. tolerance) then
                        ! Add this to a counter which keeps track of how many have converged
                        converged = converged + 1
                    end if
                end do
            end do

            ! If all points within the grid have reported that they have converged then ...
            if (converged .eq. grid_area) then
                ! Mark this victory as true
                all_converged=.true.
            end if
        end subroutine check_convergence
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to create a 2D plot of the potential
    ! across the grid, outputs image as a .pgm file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine create_pgm(Nx,Ny,array,image_name)
        implicit none
        integer, parameter :: dp = selected_real_kind(15,300)
        integer, allocatable, dimension(:,:) :: pixels
        integer :: ierr, max_greys, out_unit, Nx, Ny, i, j
        character(len=*) :: image_name

        real(kind=dp), dimension(0:Nx+1,0:Ny+1) :: array
        real(kind=dp) :: max, min

        allocate(pixels(0:Nx+1,0:Ny+1),stat=ierr)
        if (ierr.ne.0) stop "Error allocating pixels"

        ! Define values for creating the PGM
        max_greys = 255
        max = maxval(array)
        min = minval(array)

        do j = 0, Ny+1
        do i = 0, Nx+1
            ! min < T < max
            pixels(i,j) = int((array(i,j)-min)*max_greys/(max-min))
        end do
        end do

        out_unit = 10
        open(file=image_name,unit=out_unit,status='unknown',iostat=istat)
        if (istat.ne.0) stop "Error opening mpi_image.pgm"

        write(out_unit,11,iostat=istat) 'P2'                        ! PGM magic number
        if (istat.ne.0) stop "Error writing to out_unit 11"
        write(out_unit,12,iostat=istat) Nx+2, Ny+2                  ! Width & height
        if (istat.ne.0) stop "Error writing to out_unit 12"
        write(out_unit,13,iostat=istat) max_greys                   ! Max grey value
        if (istat.ne.0) stop "Error writing to out_unit 13"
        write(out_unit,fmt=*,iostat=istat) pixels(:,:)
        if (istat.ne.0) stop "Error writing to out_unit 14"

        close(unit=out_unit,iostat=istat)
        if (istat.ne.0) stop "Error closing out_unit"

        11 format(a2)
        12 format(i10,1x,i10)
        13 format(i10)

        deallocate(pixels,stat=ierr)
        if (istat.ne.0) stop "Error deallocating pixels array"
    end subroutine create_pgm

end program jacobi_serial