function reference_nonzeros (N, bandwidth) result(nonzeros)

  integer, intent(in) :: N, bandwidth
  integer :: nonzeros
  integer :: i

  nonzeros = N
  do i = 1, min(N, bandwidth)
    nonzeros = nonzeros+2*(N-i)
  enddo

end function reference_nonzeros

!> A matmul that counts the number of non-zero products.
!!
!! \f$ C \leftarrow A B \f$
!!
!! @param A Matrix A.
!! @param B Matrix B.
!! @param C Matrix C.
subroutine counting_matmul (A, B, C)

  real(kind(0d0)), intent(in) :: A(:, :), B(:, :)
  real(kind(0d0)), intent(inout) :: C(:, :)

end subroutine counting_matmul

program test

  use spammpack
  use test_utilities

#include <spamm_utility_macros.h>

  implicit none

  interface
    function reference_nonzeros (N, bandwidth) result(nonzeros)
      integer, intent(in) :: N, bandwidth
      integer :: nonzeros
    end function reference_nonzeros
  end interface

  logical, parameter :: PRINT_MATRICES = .false.
  integer, parameter :: N = 129
  integer, parameter :: BANDWIDTH = 4
  integer, parameter :: ITERATIONS = 5

  type(spamm_matrix_2nd_order), pointer :: A, B, C
  real(kind(0d0)), dimension(N, N) :: A_dense, C_dense
  real(kind(0d0)) :: reference_norm
  real(kind(0d0)), dimension(ITERATIONS) :: reference_complexity = [ 18368, 49600, 151872, 479808, 1342528 ]
  integer :: i, j

#ifndef SPAMM_COUNTERS 
  write(*, "(A)") "skipping test, library was not compiled with SPAMM_COUNTERS"
#else
  A_dense = 0
  do i = 1, N
    do j = 1, N
      if(abs(i-j) <= BANDWIDTH) then
        A_dense(i, j) = 1/dble(N)**(1/3.)
      endif
    enddo
  enddo

  if(PRINT_MATRICES) then
    call print_matrix(A_dense)
  endif

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => A
  C => spamm_zero_matrix(N, N)

  write(*, "(A)") "number non-zero = "//to_string(A%number_nonzeros)
  write(*, "(A)") "ref. non-zeros  = "//to_string(reference_nonzeros(N, BANDWIDTH))

  if(abs(A%number_nonzeros-reference_nonzeros(N, BANDWIDTH)) > 1d-10) then
    LOG_FATAL("non-zero count wrong")
    error stop
  endif

  do i = 1, ITERATIONS
    call multiply(A, B, C)

    write(*, "(A)") to_string(i)//": non-zeros       = "//to_string(C%number_nonzeros)
    write(*, "(A)") to_string(i)//": ref. non-zeros  = "//to_string(reference_nonzeros(N, 2**i*BANDWIDTH))
    write(*, "(A)") to_string(i)//": complexity      = "//to_string(C%number_operations)
    write(*, "(A)") to_string(i)//": ref. complexity = "//to_string(reference_complexity(i))
    write(*, "(A)") to_string(i)//": complexity/N^3  = "//to_string(C%number_operations/dble(N)**3)

    if(abs(C%number_nonzeros-reference_nonzeros(N, 2**i*BANDWIDTH)) > 1e-10) then
      LOG_FATAL(to_string(i)//": non-zero count wrong")
      error stop
    endif

    if(abs(C%number_operations-reference_complexity(i)) > 1e-10) then
      LOG_FATAL(to_string(i)//": complexity count wrong")
      error stop
    endif

    call copy(C, A)
    C_dense = matmul(A_dense, A_dense)
    A_dense = C_dense
  enddo

  if(PRINT_MATRICES) then
    call print_matrix(C_dense)
  endif

  reference_norm = matrix_norm(C_dense)
  write(*, "(A)") "norm_ref = "//to_string(reference_norm)
  write(*, "(A)") "norm     = "//to_string(C%norm)
  if(abs(reference_norm-C%norm) > 1d-10) then
    write(*, "(A)") "norm mismatch"
  endif

  do i = 1, size(C_dense, 1)
    do j = 1, size(C_dense, 2)
      if(abs((C_dense(i, j)-get(C, i, j))/C_dense(i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("C_reference("//to_string(i)//","//to_string(j)//") = "//to_string(C_dense(i, j)))
        LOG_FATAL("          C("//to_string(i)//","//to_string(j)//") = "//to_string(get(C, i, j)))
        error stop
      endif
    enddo
  enddo

  write(*, *) "done"
#endif

end program test
