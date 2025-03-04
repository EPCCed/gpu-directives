
! ftn -o main grid.f90 field.f90 new.f90
program new
    ! use grid_module 
    use field_module
    implicit none
    ! integer, parameter :: dp = selected_real_kind(15,300)

    type(grid_t) :: grid
    type(field_t) :: phi_old, phi_new 
    type(field_t), allocatable :: rho(:), phi1(:), phi2(:)

    integer :: i, j, k, h, iField, nFields, nIterations, nIterationsOutput
    integer :: shape(2) = [10, 10]            ! Grid shape

    real(kind=dp) :: left_box(2) = [-1.0_dp, -1.0_dp]   ! Coordinates of the bottom-left corner
    real(kind=dp) :: right_box(2) = [1.0_dp, 1.0_dp]    ! Coordinates of the top-right corner
        
    nFields = 2
    ! Simulation parameters: 
    nIterations = 100000                ! Number of time steps
    nIterationsOutput = nIterations/5   ! Number of times we output a file

    ! Create a grid using a structure (?) 
    grid = make_grid(left_box,right_box,shape)

    ! Add error checking 
    allocate(rho(nFields), phi1(nFields), phi2(nFields))

    do iField = 1, nFields
        ! Initialise grids 
        ! Allocating 0 -> 501, is that right??? 
        call rho(iField)%init(grid)
        call phi1(iField)%init(grid)
        call phi2(iField)%init(grid)

        ! Initialise fields with a gaussian 
        call init_laplacian_gaussian(10.0_dp, rho(iField), grid)
        call init_gaussian(20.0_dp, phi1(iField), grid, 0.5_dp)
        call init_gaussian(20.0_dp, phi2(iField), grid, 0.5_dp)
        
        call rho(iField)%print_field()
        print*, " "
        call phi1(iField)%print_field()
        print*, " "
        call phi2(iField)%print_field()
        print*, " "

        ! Boundary conditions -> check these are maintained by only iterating over the middle? 

        ! Print to file
        call print_to_file(rho(iField), iField, 0)
    end do

    ! Don't understand this bit with phi_new and phi_old... 


    ! do k = 1, nIterations/nIterationsOutput
    !     output_loop: do h = 1, nIterationsOutput 
    !         ! previous_grid(:,:,:) = grid(:,:,:)


    !         ! Defined in jacobi.cpp 
    !         do iField = 1, nFields
    !             do j = 2, Nx-1
    !                 do i = 2, Ny-1  
                        
    !                     ! previous_grid(i,j,iField) = 0.5*(grid(i+1,j,iField) + grid(i-1,j,iField) + grid(i,j+1,iField) + grid(i,j-1,iField) + (dx**2) * charge_density(i,j,iField))

    !                     ! grid[k - 1] = grid(i-1,j-1,iField)
    !                     nx = 
    !                     ny = 
    !                     dy = 
    !                     dx = 
    !                     aspect2 = (dy/dx)*(dy/dx)  
    !                     previous_grid(i,j,iField) = 0.5*((grid[k - 1] + grid[k + 1])/(1.0_dp + 1.0_dp/aspect2)+(grid[k - (ny+2) ] +  grid[(k + ny+2) ])/(1.0_dp + aspect2) - field_rho[k] * (dx*dx)/(1.0_dp + aspect2))
                        
    !                     ! Now update the potential value using the nearest neighbours to grid(i,j) and the charge density
    !                     !potential = 0.25_dp*(grid(i+1,j,iField) + grid(i-1,j,iField) + grid(i,j+1,iField) + grid(i,j-1,iField) + (dx**2) * charge_density(i,j,iField))
    !                     ! Update grid
    !                     !grid(i,j,iField) = grid(i,j,iField) + free_parameter*(potential - grid(i,j,iField))
    !                 end do 
    !             end do
    !         end do      
    !     end do output_loop

    !     do iField = 1, nFields
    !         call print_to_file(phi_new(iField),iField,k*nIterationsOutput)
    !     end do 
    ! end do 

    ! Add error checking 
    deallocate(rho,phi1,phi2)

contains 

subroutine init_gaussian(A, field, grid, alpha)
    implicit none
    real(kind=dp), intent(in) :: alpha, A
    type(field_t), intent(inout) :: field
    type(grid_t), intent(in) :: grid
    integer :: i, j
    real(kind=dp) :: r2
    integer :: k

    ! Loop over full grid (missing the edges for the BCs?)
    do i = 1, grid%n(1)
        do j = 1, grid%n(2) 
            ! Calculate the square of the radius r2
            r2 = grid%x(i) ** 2 + grid%y(j) ** 2
            ! Get the index in the grid
            k = grid%get_index(i,j)
            ! Assign the Gaussian value
            field%data(k) = A * exp(-alpha * r2)
        
        end do
    end do
end subroutine init_gaussian

subroutine init_laplacian_gaussian(alpha, field, grid)
    implicit none
    real(kind=dp), intent(in) :: alpha
    type(field_t), intent(inout) :: field
    type(grid_t), intent(in) :: grid
    integer :: i, j
    real(kind=dp) :: r2
    integer :: k

    ! Loop over full grid (missing the edges for the BCs?)
    do i = 1, grid%n(1)
        do j = 1, grid%n(2) 
            ! Calculate the square of the radius r2
            r2 = grid%x(i) ** 2 + grid%y(j) ** 2
            ! Get the index in the grid
            k = grid%get_index(i,j)
            ! Assign the Gaussian value
            field%data(k) = -2.0_dp * alpha*(2.0_dp - 2.0_dp * alpha * r2) * exp(-alpha * r2)
        end do
    end do
end subroutine init_laplacian_gaussian

subroutine print_to_file(field, iField, iteration)
    implicit none 
    type(field_t), intent(in) :: field
    ! type(grid_t), intent(in) :: grid
    integer :: iField, unit_num, istat, iteration
    character(len=64) :: filename_dat
    character(len=6)  :: iField_str, iteration_str
    integer(kind=8) :: nx_out, ny_out  ! Must match Python expectations
    real(8), allocatable :: field_data(:)

    write(iField_str, '(I0)') iField  
    write(iteration_str, '(I0)') iteration  

    filename_dat = "phi_" // trim(iField_str) // "_" // trim(iteration_str) // ".dat"

    unit_num = 10  
    ! open(unit=unit_num, file=filename_dat, form='unformatted', access='stream', status='unknown', iostat=istat)
    open(unit=unit_num, file=filename_dat, form='formatted', access='stream', status='unknown', iostat=istat)
    if (istat /= 0) stop "Error opening write_to_output .dat"
    
    ! Write Nx and Ny as 8-byte integers (matches Python struct.unpack('1L'))
    nx_out = field%grid%n(1)
    ny_out = field%grid%n(2)

    field_data = field%data

    write(unit_num, *, iostat=istat) nx_out
    if (istat /= 0) stop "Error writing dimensions to write_to_output .dat"
    write(unit_num, *, iostat=istat) ny_out
    if (istat /= 0) stop "Error writing dimensions to write_to_output .dat"

    ! Write grid data as raw 8-byte floats (double precision)
    ! This might be wrong
    write(unit_num, *, iostat=istat) field_data
    if (istat /= 0) stop "Error writing grid to write_to_output .dat"

    close(unit=unit_num, iostat=istat)
    if (istat /= 0) stop "Error closing write_to_output .dat"
end subroutine print_to_file

end program new 

