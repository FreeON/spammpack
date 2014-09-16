program test

  use spammpack

  implicit none

  integer, parameter :: N = 5
  integer, parameter :: LWORK = 1+6*N+2*N**2
  integer, parameter :: LIWORK = 3+5*N

  type(spamm_matrix_2nd_order), pointer :: S, X
  real(kind(0d0)), dimension(N, N) :: S_dense, X_dense
  real(kind(0d0)), dimension(N) :: eval
  real(kind(0d0)), dimension(LWORK) :: work
  integer, dimension(LIWORK) :: iwork
  integer :: info

  call random_number(S_dense)
  S_dense = S_dense+transpose(S_dense)

  S => spamm_convert_dense_to_matrix_2nd_order(S_dense)
  X => spamm_inverse_schulz(S)

#ifdef LAPACK_FOUND
  call dsyevd("V", "U", N, S_dense, N, eval, work, LWORK, iwork, LIWORK, info)
#endif

end program test
