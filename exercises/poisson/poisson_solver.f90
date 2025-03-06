program poisson_solver
    use jacobi
    implicit none

    type(grid_t) :: grid
    type(field_t), allocatable :: rho(:), phi1(:), phi2(:)

    integer :: iField, i, j, k, h                       ! Iterators
    integer :: nFields, nIterations, nIterationsOutput  ! Simulation parameters 
    integer :: Nx, Ny, index, aspect2, istat
    integer :: shape(2)                                 ! Grid shape

    real(kind=dp) :: left_box(2) = [-1.0_dp, -1.0_dp]   ! Coordinates of the bottom-left corner
    real(kind=dp) :: right_box(2) = [1.0_dp, 1.0_dp]    ! Coordinates of the top-right corner
        
    real(kind=dp) :: total_start_time, total_end_time, jacobi_start_time, jacobi_end_time, total_jacobi_time

    Nx = 500 
    Ny = 500
    shape = [Nx, Ny]

    ! Simulation parameters: 
    nFields = 1
    nIterations = 100000                ! Number of time steps
    nIterationsOutput = nIterations/5   ! Number of times we output a file

    print*, "Initialising... "

    ! Create a grid using a structure to match C++ implementation 
    grid = make_grid(left_box,right_box,shape)

    allocate(rho(nFields), phi1(nFields), phi2(nFields), stat=istat)
    if (istat.ne.0) stop "Error allocating arrays"

    do iField = 1, nFields
        ! Initialise grids 
        call rho(iField)%init(grid)
        call phi1(iField)%init(grid)
        call phi2(iField)%init(grid)

        ! Initialise fields with a gaussian 
        call init_laplacian_gaussian(10.0_dp, rho(iField), grid)
        call init_gaussian(20.0_dp, phi1(iField), grid, 0.5_dp)
        call init_gaussian(20.0_dp, phi2(iField), grid, 0.5_dp)
        
        ! Print field using: 
        ! call phi1(iField)%print_field()

        ! Boundary conditions are implict due to the method of allocation & iteration later

        call print_to_file(rho(iField), iField, 0)
    end do

    call cpu_time(total_start_time)
    print*, "Starting iterations for", nFields, "fields... "

    do k = 1, nIterations/nIterationsOutput
        do h = 1, nIterationsOutput 
            phi2 = phi1 

            call cpu_time(jacobi_start_time)
            call compute_jacobi(phi1, phi2, rho, nFields, grid)
            ! do iField = 1, nFields
            !     grid = phi1(iField)%grid

            !     do j = 2, Nx-1
            !         do i = 2, Ny-1  
            !             index = grid%get_index(i, j)
            !             aspect2 = (grid%dx(2) / grid%dx(1)) * (grid%dx(2) / grid%dx(1))
            !             phi2(iField)%data(index) = 0.5 * ((phi1(iField)%data(index - 1) + phi1(iField)%data(index + 1)) / (1.0_dp + 1.0_dp / aspect2) + (phi1(iField)%data(index - (Ny+2)) + phi1(iField)%data(index + (Ny+2))) / (1.0_dp + aspect2) - rho(iField)%data(index) * grid%dx(1) ** 2 / (1.0_dp + aspect2))
            !         end do 
            !     end do
            ! end do  
            call cpu_time(jacobi_end_time)
            total_jacobi_time = total_jacobi_time + (jacobi_end_time - jacobi_start_time)    
        end do 

        do iField = 1, nFields
            call print_to_file(phi2(iField), iField, k*nIterationsOutput)
        end do 
    end do 

    call cpu_time(total_end_time)

    print*, "Finished!"
    print*, "Took:", total_end_time-total_start_time, "s"
    print*, "On average, each jacobi took:", (total_jacobi_time/nIterations), "s"

    deallocate(rho,phi1,phi2,stat=istat)
    if (istat.ne.0) stop "Error allocating arrays"

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
    integer :: i, j, k
    real(kind=dp) :: r2

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
    integer :: iField, unit_num, istat, iteration
    character(len=64) :: filename_dat
    character(len=6) :: iField_str, iteration_str
    integer(kind=8) :: nx_out, ny_out  ! Grid sizes including ghost cells
    real(8), allocatable :: field_data(:)

    ! Convert iField and iteration to strings
    write(iField_str, '(I0)') iField
    write(iteration_str, '(I0)') iteration

    ! Construct filename for output
    filename_dat = "phi_" // trim(iField_str) // "_" // trim(iteration_str) // ".dat"

    ! Open the file for writing in binary mode
    unit_num = 10
    open(unit=unit_num,file=filename_dat,form='unformatted',access='stream',status='unknown',iostat=istat)
    if (istat.ne.0) stop "Error opening write_to_output .dat"
    
    ! Get grid dimensions and write them as 8-byte integers
    nx_out = field%grid%n(1)   ! Grid size in x-direction
    ny_out = field%grid%n(2)   ! Grid size in y-direction

    ! Write nx and ny as 8-byte integers
    write(unit_num,iostat=istat) nx_out
    if (istat.ne.0) stop "Error writing nx to write_to_output .dat"
    write(unit_num,iostat=istat) ny_out
    if (istat.ne.0) stop "Error writing ny to write_to_output .dat"

    ! Allocate field_data based on the grid size (nx_out * ny_out), including ghost cells
    allocate(field_data((nx_out + 2)*(ny_out + 2)),stat=istat)
    if (istat.ne.0) stop "Error allocating memory for field_data"

    ! Copy the data from field%data to field_data
    field_data = field%data

    ! Write field data as raw 8-byte doubles
    write(unit_num,iostat=istat) field_data
    if (istat.ne.0) stop "Error writing field data to write_to_output .dat"

    ! Close the file
    close(unit=unit_num,iostat=istat)
    if (istat.ne.0) stop "Error closing write_to_output .dat"

    deallocate(field_data,stat=istat)
    if (istat.ne.0) stop "Error allocating arrays"
end subroutine print_to_file

end program poisson_solver 

