program test

  use spammpack

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: N = 129

  type(spamm_matrix_order_2), pointer :: A => null()
  real(kind(0d0)), dimension(N, N) :: A_dense
  real(kind(0d0)) :: reference_trace
  integer :: i, j

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  reference_trace = sum(A_dense, reshape([ ((i == j, i = 1, N), j = 1, N) ], [ N, N ]))
  if(abs(reference_trace-trace(A)) > 1d-10) then
    LOG_FATAL("trace mismatch")
    error stop
  endif

end program test
