program test

  use spammpack

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: M = 13
  integer, parameter :: N = 17

  type(spamm_matrix_order_2), pointer :: A
  type(spamm_matrix_order_1), pointer :: V
  real(kind(0d0)), dimension(M, N) :: A_dense
  integer :: i
  real(kind(0d0)) :: Aij, Vi

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  call spamm_allocate_matrix_order_1(M, V)

  call copy(A, N, V)

  do i = 1, M
    Aij = get(A, i, N)
    Vi = get(V, i)
    if(abs(Aij-Vi) /= 0) then
      LOG_FATAL("A("//to_string(i)//","//to_string(N)//") mismatch")
      LOG_FATAL("           found "//to_string(Vi))
      LOG_FATAL("should have been "//to_string(Aij))
      error stop
    endif
  enddo

  call delete(A)
  call delete(V)

end program test
