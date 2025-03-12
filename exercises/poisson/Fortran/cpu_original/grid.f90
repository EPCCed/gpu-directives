module grid_module
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)

    type :: grid_t
        integer :: n(2)      ! Grid dimensions
        real(kind=dp) :: dx(2)     ! Grid spacing in x and y
        real(kind=dp) :: start(2)  ! Start coordinates
        real(kind=dp) :: end(2)    ! End coordinates
    contains
        procedure :: get_index
        procedure :: x
        procedure :: y
        procedure :: size
    end type grid_t

contains

    ! Returns the index in memory associated with cell (i, j)
    pure function get_index(this, i, j) result(index)
        class(grid_t), intent(in) :: this
        integer, intent(in) :: i, j
        integer :: index
        index = (j + 1) + (i + 1) * (this%n(2) + 2)
    end function get_index
    
    ! Computes x-coordinate at index i
    pure function x(this, i) result(x_coord)
        class(grid_t), intent(in) :: this
        integer, intent(in) :: i
        real(8) :: x_coord
        x_coord = this%start(1) + this%dx(1) * i
    end function x

    ! Computes y-coordinate at index j
    pure function y(this, j) result(y_coord)
        class(grid_t), intent(in) :: this
        integer, intent(in) :: j
        real(8) :: y_coord
        y_coord = this%start(2) + this%dx(2) * j
    end function y

    ! Computes the total number of grid cells including ghost cells
    pure function size(this) result(grid_size)
        class(grid_t), intent(in) :: this
        integer :: grid_size
        grid_size = (this%n(1) + 2) * (this%n(2) + 2)
    end function size

    ! Creates and initializes a grid
    function make_grid(start, end, n) result(grid)
        real(kind=dp), intent(in) :: start(2), end(2)
        integer, intent(in) :: n(2)
        type(grid_t) :: grid

        grid%n = n
        grid%start = start
        grid%end = end
        grid%dx(1) = (end(1) - start(1)) / n(1)
        grid%dx(2) = (end(2) - start(2)) / n(2)
    end function make_grid

end module grid_module
