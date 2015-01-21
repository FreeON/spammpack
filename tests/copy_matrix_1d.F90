program copy_matrix

  use spammpack

  implicit none

  type(spamm_tree_1d) :: A, B

#ifdef HAVE_CONSTRUCTOR
  A = spamm_tree_1d(10)
#else
  A = new_tree_1d(10)
#endif
  B = A

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

end program copy_matrix
