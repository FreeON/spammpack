program multiply

  use spammpack

  type(spamm_matrix_2nd_order), pointer :: A, B, C
  real(spamm_kind), dimension(5, 5) :: A_dense, C_dense

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  C_dense = 0

end program multiply
