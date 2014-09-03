program test

#include "spamm_utility_macros.h"

  use spammpack

  implicit none

  integer, parameter :: M = 102
  integer, parameter :: N = 78

  type(spamm_matrix_2nd_order), pointer :: A => null()
  type(spamm_matrix_order_1), pointer :: B => null()
  type(spamm_matrix_order_1), pointer :: C => null()

  real(kind(0d0)), dimension(M, N) :: A_dense
  real(kind(0d0)), dimension(N) :: B_dense
  real(kind(0d0)), dimension(M) :: C_dense

  integer :: i

  call random_number(A_dense)
  call random_number(B_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => spamm_convert_dense_to_order_1(B_dense)

  call multiply(A, B, C)

  C_dense = matmul(A_dense, B_dense)

  do i = 1, M
    if(abs(get(C, i)-C_dense(i)) > 1d-10) then
      LOG_FATAL("vector element mismatch")
      error stop
    endif
  enddo

end program test
