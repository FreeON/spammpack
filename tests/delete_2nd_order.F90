program test

  use spammpack

  implicit none

  integer, parameter :: M = 5
  integer, parameter :: N = 6

  type(spamm_matrix_order_2), pointer :: A
  real(kind(0d0)), dimension(M, N) :: A_dense

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  call delete(A)

end program test
