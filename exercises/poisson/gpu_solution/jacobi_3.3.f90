module jacobi
    use field_module
    use omp_lib
    implicit none 

    contains 
    subroutine compute_jacobi(phi1, phi2, rho, nFields, grid)
        implicit none 
        integer :: iField, nFields, i, j, Nx, Ny 
        integer :: index = 0, aspect2 = 0 

        ! Argument declarations
        type(field_t), dimension(:), intent(inout) :: phi1, phi2, rho
        type(grid_t) :: grid
        Nx = grid%n(1)
        Ny = grid%n(2)

        ! Parallelise across multiple CPU threads. 
        ! Each thread will handle different iterations of the loop.
        !$omp parallel do 
        do iField = 1, nFields
            grid = phi1(iField)%grid
            ! Offload this bit to the GPU: 
            !$omp target teams distribute parallel do collapse(2) 
            do j = 2, Nx-1
                do i = 2, Ny-1  
                    ! Inlined get index
                    index = (j + 1) + (i + 1) * (grid%n(2) + 2)                    
                    aspect2 = (grid%dx(2) / grid%dx(1)) * (grid%dx(2) / grid%dx(1))
                    phi2(iField)%data(index) = 0.5 * ((phi1(iField)%data(index - 1) + phi1(iField)%data(index + 1)) / (1.0_dp + 1.0_dp / aspect2) + (phi1(iField)%data(index - (Ny+2)) + phi1(iField)%data(index + (Ny+2))) / (1.0_dp + aspect2) - rho(iField)%data(index) * grid%dx(1) ** 2 / (1.0_dp + aspect2))
                end do 
            end do
            !$omp end target teams distribute parallel do
        end do  
        !$omp end parallel do 

    end subroutine compute_jacobi 

end module jacobi