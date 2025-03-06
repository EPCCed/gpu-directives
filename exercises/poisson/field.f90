module field_module
    use grid_module
    implicit none

    type :: field_t
        real(8), allocatable :: data(:)   ! Use allocatable instead of pointer
        type(grid_t) :: grid              ! Direct storage instead of a pointer
    contains
        procedure :: init
        procedure :: get_data
        procedure :: get_grid
        procedure :: print_field
    end type field_t

contains

    ! Initializes the field based on the given grid
    subroutine init(this, grid_)
        class(field_t), intent(inout) :: this
        type(grid_t), intent(in) :: grid_

        this%grid = grid_  ! Copy grid 
        if (allocated(this%data)) deallocate(this%data)
        allocate(this%data(grid_%size()))  ! Allocate memory for data
        this%data = 0.0d0  ! Initialize to zero
    end subroutine init

    ! Returns field data (direct access)
    function get_data(this) result(data_out)
        class(field_t), intent(in) :: this
        real(8), allocatable :: data_out(:)
        data_out = this%data  
    end function get_data

    ! Returns associated grid (direct access)
    function get_grid(this) result(grid_out)
        class(field_t), intent(in) :: this
        type(grid_t) :: grid_out
        grid_out = this%grid  ! Copy grid (no pointers involved)
    end function get_grid

    subroutine print_field(this)
        class(field_t), intent(in) :: this
        integer :: i, j, k
        do i = 0, this%grid%n(1) + 1
            ! Print each row's values in one line
            do j = 0, this%grid%n(2) + 1
                k = this%grid%get_index(i, j)
                write(*, "(F8.3)", advance="no") this%data(k)
                ! Optionally, you can add a space between values
                if (j < this%grid%n(2) + 1) then
                    write(*, "(A)", advance="no") " "  ! Add space between columns
                end if
            end do
            print *  ! Move to the next line after printing all values in the row
        end do
    end subroutine print_field
end module field_module
