program copy_matrix

  use spammpack

  implicit none

  type(tree_1d), pointer :: A, B

  A => new_tree_1d(10)
  B => A

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

end program copy_matrix
