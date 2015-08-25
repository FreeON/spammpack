program test

  use spammpack
  implicit none

  integer, parameter :: N = 100
  double precision, parameter :: alpha = 1.2
  double precision, parameter :: beta = 0.8

  type(spamm_tree_2d_symm), pointer :: a
  type(spamm_tree_2d_symm), pointer :: b

  double precision :: a_dense(N, N)
  double precision :: b_dense(N, N)

  call random_number(a_dense)
  call random_number(b_dense)
  a => spamm_convert_from_dense(a_dense)
  b => spamm_convert_from_dense(b_dense)

  a_dense = alpha*a_dense+beta*b_dense
  a => spamm_tree_2d_symm_plus_tree_2d_symm(a, b, alpha, beta, a)

  call spamm_convert_tree_2d_symm_to_dense(a, b_dense)
  if(maxval(a_dense-b_dense) > 1d-10) then
     write(*, *) "Value mismatch"
     error stop
  end if
  write(*, *) "matrices match"

  call spamm_destruct_tree_2d_symm_recur(a)
  call SpAMM_destruct_tree_2d_symm_recur(b)

end program test
