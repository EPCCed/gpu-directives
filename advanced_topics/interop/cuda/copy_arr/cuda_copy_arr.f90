Program main

  use iso_c_binding

  implicit none

  integer, parameter :: fp_kind = kind(0.0)
  type(c_ptr) :: cptr_b, cptr_c
  real(fp_kind), pointer, dimension(:) :: fptr_b => null(), fptr_c => null()
  integer :: n = 5, i, res
  real, allocatable, dimension(:), target :: b_h, c_h
  integer * 8 cudaMemcpyDeviceToHost, cudaMemcpyHostToDevice
  parameter(cudaMemcpyHostToDevice=1)
  parameter(cudaMemcpyDeviceToHost=2)

  interface
    integer(C_INT) function cudaMalloc(buffer, size) bind(C, name="cudaMalloc")
      use iso_c_binding
      implicit none
      type(C_PTR) :: buffer
      integer(C_SIZE_T), value :: size
    end function cudaMalloc
    integer function cudaMemcpy(dst, src, count, kdir) bind(C, name="cudaMemcpy")
      ! note: cudaMemcpyHostToDevice = 1
      ! note: cudaMemcpyDeviceToHost = 2
      use iso_c_binding
      implicit none
      type(C_PTR), value :: dst, src
      integer(c_size_t), value :: count, kdir
    end function
  end interface

  allocate (b_h(n), c_h(n))

  do i = 1, n
    b_h(i) = i
  end do

  print *, b_h

  res = cudaMalloc(cptr_b, n * sizeof(fp_kind))
  res = cudaMalloc(cptr_c, n * sizeof(fp_kind))
  res = cudaMemcpy(cptr_b, c_loc(b_h), n * sizeof(fp_kind), cudaMemcpyHostToDevice)

  call c_f_pointer(cptr_b, fptr_b, [n])
  call c_f_pointer(cptr_c, fptr_c, [n])

!$omp target teams distribute parallel do simd is_device_ptr(fptr_b, fptr_c)
  do i = 1, n
    fptr_c(i) = fptr_b(i)
  end do

  res = cudaMemcpy(c_loc(c_h), cptr_c, n * sizeof(fp_kind), cudaMemcpyDeviceToHost)

  print *, c_h

end program main
