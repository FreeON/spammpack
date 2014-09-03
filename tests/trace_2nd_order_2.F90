program test

#include "spamm_utility_macros.h"

  use spammpack

  implicit none

  integer, parameter :: N = 129

  type(spamm_matrix_2nd_order), pointer :: A => null()
  type(spamm_matrix_2nd_order), pointer :: B => null()

  real(kind(0d0)), dimension(N, N) :: A_dense, B_dense, C_dense

  real(kind(0d0)) :: reference_trace
  integer :: i, j

  call random_number(A_dense)
  call random_number(B_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => spamm_convert_dense_to_matrix_2nd_order(B_dense)

  C_dense = matmul(A_dense, B_dense)

  reference_trace = sum(C_dense, reshape([ ((i == j, i = 1, N), j = 1, N) ], [ N, N ]))

  if(abs(reference_trace-trace(A, B)) > 1d-10) then
    LOG_FATAL("trace mismatch")
    error stop
  endif

end program test
