! Serial version
program jacobi_serial
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    integer :: i, j, k, h, istat
    integer :: Nx, Ny, Lx, Ly, dx, dy, x1, x2, grid_area
    integer :: nIterations, nIterationsOutput, nFields, iField
    
    real(kind=dp) :: A, alpha, r2
    real(kind=dp) :: potential, free_parameter, tolerance
    real(kind=dp), allocatable, dimension(:,:,:) :: grid, previous_grid, charge_density
    
    logical :: all_converged
    real(kind=dp) :: start_time, finish_time, start_jacobi_time, finish_jacobi_time, total_jacobi_time  

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
    
    nFields = 2                     ! Number of field equations to solve

    ! Allocate a Nx * Ny grid for each of the field equations
    allocate(grid(0:Nx+1,0:Ny+1,nFields),stat=istat)
    if (istat.ne.0) stop "Error allocating grid array"
    allocate(previous_grid(0:Nx+1,0:Ny+1,nFields),stat=istat)
    if (istat.ne.0) stop "Error allocating previous_grid array"
    allocate(charge_density(0:Nx+1,0:Ny+1,nFields),stat=istat)
    if (istat.ne.0) stop "Error allocating charge_density array"

    ! Simulation parameters: 
    nIterations = 100000                ! Number of time steps
    nIterationsOutput = nIterations/5   ! Number of times we output a file
    free_parameter = 0.8_dp             ! Overelaxation factor (must be less than one for overrelaxation)
    tolerance = 1.0E-9_dp               ! Tolerance to accept convergence
    start_time = 0.0_dp 
    finish_time = 0.0_dp
    
    ! Gaussian distribution 
    A = 0.1_dp                      ! Max value of potential
    alpha = 0.01_dp                 ! Width of gaussian 

    print*, "Initialising..."

    do iField = 1, nFields
        grid(:,:,iField) = 0.0_dp 
        previous_grid(:,:,iField) = 0.0_dp 
        charge_density(:,:,iField) = 0.0_dp 

        do i = 1, Nx
            do j = 1, Ny
                r2 = (real(i-x1, dp)**2 + real(j-x2, dp)**2)
                ! Initialise the field with either 0 or gaussian -> This helps for a smooth convergence.
                grid(i,j,iField) = A * exp(-alpha * r2)
                ! grid(i,j,2) = 0.2_dp * exp(-alpha * r2)
                ! Initialise the charge density with a uniform or gaussian distribution. 
                charge_density(i,j,iField) = A * (2.0_dp * alpha - 4.0_dp * alpha**2 * r2) * exp(-alpha * r2)
                ! charge_density(i,j,iField) = 3.0_dp 
            end do
        end do

        ! Output the initial conditions
        call write_to_output(0,Nx,Ny,grid,iField)
    end do 

    call cpu_time(start_time)

    print*, "Starting calculation with", nFields, "fields ..."
               
    ! Because we are outputting every nIterationsOutput, we need to do nIterations/nIterationsOutput chunks to complete
    ! the total number of iterations
    do k = 1, nIterations/nIterationsOutput
        output_loop: do h = 1, nIterationsOutput 
            ! Store the previous grid
            previous_grid(:,:,:) = grid(:,:,:)

            ! Jacobi: 
            call cpu_time(start_jacobi_time)

            ! The grid is initialised with 0 values. By iterating between 1 and Nx/Ny - 1, we are maintaining a boundary of 0, 
            ! therefore maintaining Dirichlet boundary conditions. 

            !!!!!!!!!!!!!!!!!!!!!!! omp would go round here: 

            do iField = 1, nFields
                do j = 2, Nx-1
                    do i = 2, Ny-1  
                        ! Now update the potential value using the nearest neighbours to grid(i,j) and the charge density
                        potential = 0.25_dp*(grid(i+1,j,iField) + grid(i-1,j,iField) + grid(i,j+1,iField) + grid(i,j-1,iField) + (dx**2) * charge_density(i,j,iField))
                        ! Update grid
                        grid(i,j,iField) = grid(i,j,iField) + free_parameter*(potential - grid(i,j,iField))
                    end do 
                end do
            end do      
            !!!!!!!!!!!!!!!!!!!!!!!
            
            call cpu_time(finish_jacobi_time)
            total_jacobi_time = total_jacobi_time + (finish_jacobi_time - start_jacobi_time)
        end do output_loop

        do iField = 1, nFields
            call write_to_output(k*nIterationsOutput,Nx,Ny,grid,iField)
        end do 
    end do 

    call cpu_time(finish_time)

    print*, "Finalising..."
    print*, " "
    print*, "Average time spent in jacobi:", (total_jacobi_time/nIterations)/1000.0_dp, "ms/it"
    print*, "Total time to compute", nFields, "fields:", finish_time-start_time, "s"
    print*, " "

    deallocate(grid,stat=istat)
    if (istat.ne.0) stop 'Error deallocating grid array'
    deallocate(previous_grid,stat=istat)
    if (istat.ne.0) stop 'Error deallocating previous_grid array'
    deallocate(charge_density,stat=istat)
    if (istat.ne.0) stop 'Error deallocating previous_grid array'
  
contains 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the grid has converged
    ! by comparing the previous value point at the current grid point
    ! Usage: 
    !       call check_convergence 
    !       if (all_converged) print*, "Field", iField, "finished in", k*h, "iterations"
    !       if (all_converged) call write_to_output(k*h,Nx,Ny,grid,iField)
    !       if (all_converged) exit output_loop 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine check_convergence
    !     implicit none
    !     integer,parameter :: dp = selected_real_kind(15,300)
    !     integer :: converged 
        
    !     converged = 0
    !     all_converged = .false.
    !     ! This is an expensive test
    !     do i = 1, Nx
    !         do j = 1, Ny
    !             ! Check that each point value is within some tolerance of the previous value of that point
    !             if (abs(grid(i,j,iField)-previous_grid(i,j,iField)) .lt. tolerance) then
    !                 ! Add this to a counter which keeps track of how many have converged
    !                 converged = converged + 1
    !             end if
    !         end do
    !     end do

    !     ! If all points within the grid have reported that they have converged then ...
    !     if (converged .eq. grid_area) then
    !         ! Mark this victory as true
    !         all_converged=.true.
    !     end if
    ! end subroutine check_convergence

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to write to output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_to_output(k,Nx,Ny,grid,iField)
        implicit none
        integer, intent(in) :: k, Nx, Ny, iField
        real(kind=8), intent(in) :: grid(0:Nx+1, 0:Ny+1, iField)  
        character(len=64) :: filename_dat, filename_pgm 
        character(len=6)  :: k_str, iField_str 
        integer :: istat, unit_num
        integer(kind=8) :: nx_out, ny_out  ! Must match Python expectations

        ! Convert integer k to string
        write(k_str, '(I0)') k  
        write(iField_str, '(I0)') iField  

        ! Construct the output filename
        filename_dat = "phi_" // trim(k_str) // "_field_" // trim(iField_str) // ".dat"
        filename_pgm = "phi_" // trim(k_str) // "_field_" // trim(iField_str) // ".pgm"
    
        ! Visualise output 
        call create_pgm(Nx,Ny,grid,filename_pgm,iField)

        ! Open binary file with stream access (avoids extra record markers)
        unit_num = 10  
        open(unit=unit_num, file=filename_dat, form='unformatted', access='stream', status='unknown', iostat=istat)
        if (istat /= 0) stop "Error opening write_to_output .dat"

        ! Write Nx and Ny as 8-byte integers (matches Python struct.unpack('1L'))
        nx_out = Nx
        ny_out = Ny
        write(unit_num, iostat=istat) nx_out, ny_out
        if (istat /= 0) stop "Error writing dimensions to write_to_output .dat"

        ! Write grid data as raw 8-byte floats (double precision)
        write(unit_num, iostat=istat) grid(:,:,iField)
        if (istat /= 0) stop "Error writing grid to write_to_output .dat"

        ! Close the file
        close(unit=unit_num, iostat=istat)
        if (istat /= 0) stop "Error closing write_to_output .dat"

    end subroutine write_to_output
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to create a 2D plot of the potential
    ! across the grid, outputs image as a .pgm file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine create_pgm(Nx,Ny,array,image_name,iField)
        implicit none
        integer, parameter :: dp = selected_real_kind(15,300)
        integer, allocatable, dimension(:,:) :: pixels
        integer :: ierr, max_greys, out_unit, Nx, Ny, i, j, iField
        character(len=*) :: image_name

        real(kind=dp), dimension(0:Nx+1,0:Ny+1,iField) :: array
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
                pixels(i,j) = int((array(i,j,iField)-min)*max_greys/(max-min))
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
        ! Write pixel data in a safe, formatted way
        do j = 0, Ny+1
            write(out_unit,'(1000I5)',iostat=istat) pixels(:, j)
            if (istat.ne.0) stop "Error writing to out_unit 13"
        end do
            
        close(unit=out_unit,iostat=istat)
        if (istat.ne.0) stop "Error closing out_unit"

        11 format(a2)
        12 format(i10,1x,i10)
        13 format(i10)

        deallocate(pixels,stat=ierr)
        if (istat.ne.0) stop "Error deallocating pixels array"
    end subroutine create_pgm

end program jacobi_serial