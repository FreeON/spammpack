program test

  use spammpack

  implicit none

  type(spamm_matrix_2nd_order), pointer :: A
  real(kind(0d0)), dimension(10, 10) :: A_dense

  call random_number(A_dense)
  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

end program test
