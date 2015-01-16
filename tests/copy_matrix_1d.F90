program copy_matrix

  use spammpack

  implicit none

  type(spamm_matrix_1d) :: A, B

  A = spamm_matrix_1d(10)
  B = A

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

end program copy_matrix
