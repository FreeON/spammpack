program test

  use spammpack

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: M = 5
  integer, parameter :: N = 6

  type(spamm_matrix_order_2), pointer :: A, B
  real(kind(0d0)), dimension(M, N) :: A_dense
  integer :: i, j
  real(kind(0d0)) :: Aij, Bij

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => null()

  call copy(A, B)

  do i = 1, M
    do j = 1, N
      Aij = get(A, i, j)
      Bij = get(B, i, j)
      if(abs(Aij-Bij) /= 0) then
        LOG_FATAL("A("//to_string(i)//","//to_string(j)//") mismatch")
        LOG_FATAL("           found "//to_string(Bij))
        LOG_FATAL("should have been "//to_string(Aij))
        error stop
      endif
    enddo
  enddo

  call delete(A)
  call delete(B)

end program test
