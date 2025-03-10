Program main

  implicit none

  real*8, allocatable, dimension(:, :) :: A, B, C
  integer :: N = 4, M = 2, K = 3

  interface
    subroutine cublasDgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) &
      bind(c, name='cublasDgemm')
      use iso_c_binding
      integer(c_int), value :: m, n, k, lda, ldb, ldc
      real(c_double), intent(in), device :: a(lda, *), b(ldb, *)
      real(c_double), intent(inout), device :: c(ldc, *)
      real(c_double), value :: alpha, beta
      character(kind=c_char), value :: transa, transb
    end subroutine
    subroutine cublasInit() bind(c, name="cublasInit")
    end subroutine
    subroutine cublasShutdown() bind(c, name="cublasShutdown")
    end subroutine
  end interface

  allocate (A(M, K))
  allocate (B(K, N))
  allocate (C(M, N))

  A = 0
  B = 0
  C = 0

  call random_number(A)
  call random_number(B)

  call cublasInit

  !$omp target data map(to:A(1:M,1:K),B(1:K,1:N)) map(tofrom:C(1:M,1:N))
  !$omp target data use_device_addr(A(1:M,1:K), B(1:K,1:N), C(1:M,1:N))

  call cublasDgemm('N', 'N', M, N, K, 1.0d0, A, M, B, K, 0.0d0, C, M)

  !$omp end target data
  !$omp end target data

  call cublasShutdown

  print *, "A"
  call print_matrix(A, M, K)
  print *, ""
  print *, "B"
  call print_matrix(B, K, N)
  print *, ""
  print *, "C"
  call print_matrix(C, M, N)

end program main

subroutine print_matrix(matrix, rows, cols)
  implicit none
  integer, intent(in) :: rows, cols
  real*8, intent(in) :: matrix(rows, cols)
  integer :: i, j

  do i = 1, rows
    do j = 1, cols
      write (*, '(F8.6)', advance='no') matrix(i, j)
      if (j < cols) then
        write (*, '(A)', advance='no') "  "
      end if
    end do
    write (*, *)
  end do
end subroutine print_matrix
