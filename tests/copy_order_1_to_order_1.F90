program test

  use spammpack

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: N = 6

  type(spamm_matrix_order_1), pointer :: A, B
  real(kind(0d0)), dimension(N) :: A_dense
  integer :: i
  real(kind(0d0)) :: Ai, Bi

  call random_number(A_dense)

  A => spamm_convert_dense_to_order_1(A_dense)
  B => null()

  call copy(A, B)

  do i = 1, N
    Ai = get(A, i)
    Bi = get(B, i)
    if(abs(Ai-Bi) /= 0) then
      LOG_FATAL("A("//to_string(i)//") mismatch")
      LOG_FATAL("           found "//to_string(Bi))
      LOG_FATAL("should have been "//to_string(Ai))
      error stop
    endif
  enddo

  call delete(A)
  call delete(B)

end program test
