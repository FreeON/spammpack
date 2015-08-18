program testprogram

  use spammpack

  implicit none

  integer, parameter :: N = 10
  type(tree_2d_symmetric), pointer :: A

  A => identity_tree_2d_symmetric(N)

  write(*, "(A)") "A: "//trim(A%to_string())

  call delete_tree_2d_symmetric(A)

end program testprogram
