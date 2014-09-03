program test

  use spammpack

  implicit none

  integer, parameter :: M = 102
  integer, parameter :: N = 78

  type(spamm_matrix_2nd_order), pointer :: A => null()
  type(spamm_matrix_order_1), pointer :: V => null()

  real(kind(0d0)), dimension(M, N) :: A_dense
  real(kind(0d0)), dimension(M) :: V_dense

  call random_number(A_dense)
  call random_number(V_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  V => spamm_convert_dense_to_order_1(V_dense)

  error stop

end program test
