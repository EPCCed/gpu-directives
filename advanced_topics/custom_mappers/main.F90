
module my_structure_mod 
    use iso_fortran_env, only:  dp=>real64
    implicit none 

    type :: my_type

        real(kind=dp) :: value
        real(kind=dp) , allocatable :: c(:)
        integer :: n

        contains 
            procedure :: add_gpu
    end type

    contains
        subroutine add_gpu( this )
            class(my_type), intent(inout) :: this
            integer :: i 

!$omp target teams distribute parallel do
            do i=1,this%n 
                this%c(i)=this%c(i) + this%value
            end do
        end subroutine 

end module 


program main 
    use iso_fortran_env, only:  dp=>real64
    use my_structure_mod
    implicit none
    real(kind= dp) :: expected
    real(kind=dp) , parameter :: tol = 1e-6
    integer :: i 

    !$omp declare mapper(my_type :: d ) map(d,d%value,d%n,d%c(1:d%n))


    type(my_type) :: a 

    a%value = 2
    a%n = 10
    allocate( a%c(1:a%n) )
    
    do i=1,a%n
        a%c(i)=i 
    end do
!$omp target data map(tofrom:a)
        call a%add_gpu();
!$omp end target data


    do i=1,a%n
        expected = ( i + a%value );
        if  ( abs(a%c(i) - expected) > tol ) then

            print *, "Got ", a%c(i) ,  "instead of" , expected
            call exit(1)
        endif
    end do


    print * ,"Success!"
    contains




end program main