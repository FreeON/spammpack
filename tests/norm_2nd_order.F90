program test

  use spammpack

  implicit none

  integer, parameter :: M = 5
  integer, parameter :: N = 11

  type(spamm_matrix_2nd_order), pointer :: A
  real(kind(0d0)), dimension(M, N) :: A_dense
  real(kind(0d0)) :: norm_reference, norm_A

  call random_number(A_dense)
  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  norm_reference = sqrt(sum(A_dense**2))
  norm_A = sqrt(norm(A))

  if(abs(norm_A-norm_reference) > 1d-12) then
    write(*, *) "norm mismatch"
    write(*, *) "norm_reference = ", norm_reference
    write(*, *) "       norm(A) = ", norm_A
    error stop
  endif

  call delete(A)

end program test
