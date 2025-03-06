module jacobi
    use field_module
    implicit none 

    contains 
    subroutine compute_jacobi(phi1, phi2, rho, nFields, grid)
        implicit none 
        integer :: iField, nFields, i, j, Nx, Ny 
        integer :: index, aspect2

        ! Argument declarations
        type(field_t), dimension(:), intent(inout) :: phi1, phi2, rho
        type(grid_t) :: grid
        Nx = grid%n(1)
        Ny = grid%n(2)

        do iField = 1, nFields
            grid = phi1(iField)%grid

            do j = 2, Nx-1
                do i = 2, Ny-1  
                    index = grid%get_index(i, j)
                    aspect2 = (grid%dx(2) / grid%dx(1)) * (grid%dx(2) / grid%dx(1))
                    phi2(iField)%data(index) = 0.5 * ((phi1(iField)%data(index - 1) + phi1(iField)%data(index + 1)) / (1.0_dp + 1.0_dp / aspect2) + (phi1(iField)%data(index - (Ny+2)) + phi1(iField)%data(index + (Ny+2))) / (1.0_dp + aspect2) - rho(iField)%data(index) * grid%dx(1) ** 2 / (1.0_dp + aspect2))
                end do 
            end do
        end do  

    end subroutine compute_jacobi 

end module jacobi