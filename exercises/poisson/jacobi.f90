! Goal is: Find the electrostatic potential inside a 2D box with earthed conducting walls
! a charge applied across it? by solving Poisson's equation with appropriate boundary conditions. 

! Using the SOR method: https://people.eecs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html
! Term1_HPC/T1Exam_Q/FinalAssignment_PartB

! TO DO: 
! * Check if I need the h spacing parameter... Apparently it cancelled out in my assignment, but why? 
!       Using a grid with grid spacing (h). U(i,j) is the approximate solution at x = i*h and y = j*h. 
!       We solve U(i,j) using it's neighbours. 
!       Boundaries are 0. 
!       Charge density = b(i,j) = -f(i*h,j*h)*h^2, 
! * Check if this is working as expected? 
! * Apply a gaussian to the charge_density? Currently just a random fixed number 
! * Experiment with free_parameter to optimize for convergence (values around 0.8-1.0).

! Serial version
program jacobi_serial
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    integer :: nIterations, i, j, k, Nx, Ny, istat, converged, grid_area
    real(kind=dp), allocatable, dimension(:,:) :: grid, previous_grid, charge_density
    real(kind=dp) :: potential, free_parameter, tolerance, A, alpha
    logical :: all_converged

    Nx = 100                    ! Dimension in x
    Ny = 100                    ! Dimension in y 
    nIterations = 100000        ! Number of time steps
    free_parameter = 0.8_dp     ! Overelaxation factor (must be less than one for overrelaxation)
    tolerance = 1.0E-6_dp       ! Tolerance to accept convergence
    grid_area = Nx*Ny           
    ! h = 1.0_dp/real(scale,kind=dp)
    
    A = 1.0_dp                  ! Amplitude of gaussian
    alpha = 0.001_dp            ! Width of gaussian 


    allocate(grid(1:Nx,1:Ny),stat=istat)
    if (istat.ne.0) stop "Error allocating grid array"
    allocate(previous_grid(1:Nx,1:Ny),stat=istat)
    if (istat.ne.0) stop "Error allocating previous_grid array"
    allocate(charge_density(1:Nx,1:Ny),stat=istat)
    if (istat.ne.0) stop "Error allocating charge_density array"

    grid(:,:) = 0.0_dp 
    previous_grid(:,:) = 0.0_dp 
    charge_density(:,:) = 1.0_dp 

    ! Initialise the field? 
    do i = 1, Nx
        do j = 1, Ny
            grid(i,j) = A * exp(-alpha * (real(i, dp)**2 + real(j, dp)**2))
        end do
    end do
    
    ! Dirichlet boundary conditions = setting the outer edges to a constant value. 
    ! i.e. grid (i,1), grid(i,Ny), grid(1,j), grid(Nx,j) = 0
    do j = 1, Ny
        grid(1,j) = 0.0_dp
        grid(Nx,j) = 0.0_dp 
    end do 
    do i = 1, Nx
        grid(i,1) = 0.0_dp 
        grid(i,Ny) = 0.0_dp 
    end do 

    ! Uncomment to print the grid
    ! do j = 1, Ny 
    !     print*, grid(:,j)
    ! end do 
    
    convergence_loop: do k = 1, nIterations 
        ! Store the previous grid
        previous_grid(:,:) = grid(:,:)

        ! Reset parameters for each iteration
        converged = 0

        ! Maintain boundary conditions
        do j = 2, Nx-1
            do i = 2, Ny-1  
                ! Now update the potential value using the nearest neighbours to grid(i,j) and the charge density
                potential = 0.25_dp*(grid(i+1,j) + grid(i-1,j) + grid(i,j+1) + grid(i,j-1) + charge_density(i,j))
                ! Update grid
                grid(i,j) = grid(i,j) + free_parameter*(potential - grid(i,j))
            end do 
        end do     
    
        call check_convergence 
        if (all_converged) print*, "Converged"
        if (all_converged) exit convergence_loop
    end do convergence_loop

    deallocate(grid,stat=istat)
    if (istat.ne.0) stop 'Error deallocating grid array'
    deallocate(previous_grid,stat=istat)
    if (istat.ne.0) stop 'Error deallocating previous_grid array'
  
contains 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to check if the grid has converged
  ! by comparing the previous value point at the current grid point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_convergence
        implicit none
        integer,parameter :: dp = selected_real_kind(15,300)

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

        print*, "Converged = ", converged, "  Target = ", grid_area

        ! If all points within the grid have reported that they have converged then ...
        if (converged .eq. grid_area) then
            ! Mark this victory as true
            all_converged=.true.
        end if
    end subroutine check_convergence
end program jacobi_serial