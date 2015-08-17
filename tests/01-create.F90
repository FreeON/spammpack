program test

  use spamm_m

  implicit none

  integer, parameter :: N = 100
  type(spamm_tree_2d), pointer :: A

  A => new_tree_2d(N)

end program test
