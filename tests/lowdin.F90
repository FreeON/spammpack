program test

  use spammpack
  use test_utilities

  implicit none

#ifdef LAPACK_FOUND
  interface
    subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LIWORK, LWORK, N
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
    end subroutine dsyevd
  end interface
#endif

  integer, parameter :: N = 253
  integer, parameter :: LWORK = 1+6*N+2*N**2
  integer, parameter :: LIWORK = 3+5*N

  type(spamm_matrix_order_2), pointer :: S, Y, Z
  real(kind(0d0)), dimension(N, N) :: S_dense, Y_dense, Z_dense
  real(kind(0d0)), dimension(N) :: eval
  real(kind(0d0)), dimension(N, N) :: evec
  real(kind(0d0)), dimension(LWORK) :: work
  integer, dimension(LIWORK) :: iwork
  integer :: info
  integer :: i, j, k
  real(kind(0d0)) :: abs_diff, max_diff, Y_norm, Z_norm

  call random_number(eval)
  call random_number(evec)

  eval = eval/2+0.4

  do k = 1, N
    evec(k, :) = evec(k, :)/sqrt(sum(evec(k, :)*evec(k,:)))
  enddo

  !call print_matrix_python(reshape(eval, [ 1, size(eval) ]), "lambda")
  !call print_matrix_python(evec, "basis")

  S_dense = 0
  do k = 1, N
    forall(i = 1:N, j = 1:N)
      S_dense(i, j) = S_dense(i, j)+eval(k)*evec(k, i)*evec(k, j)
    end forall
  enddo

  !call print_matrix_python(S_dense, "S")

  S => spamm_convert_dense_to_matrix_2nd_order(S_dense)
  Y => null()
  Z => null()
  call spamm_inverse_sqrt_schulz(S, Y, Z)

#ifdef LAPACK_FOUND
  call dsyevd("V", "U", N, S_dense, N, eval, work, LWORK, iwork, LIWORK, info)

  write(*, "(A)") "condition number: "//to_string(eval(N)/eval(1))

  !call print_matrix_python(reshape(eval, [ 1, size(eval) ]), "eigenvalue")
  !call print_matrix_python(S_dense, "eigenvector")

  Y_dense = 0
  Z_dense = 0
  do k = 1, N
    forall(i = 1:N, j = 1:N)
      Y_dense(i, j) = Y_dense(i, j)+sqrt(eval(k))*S_dense(i, k)*S_dense(j, k)
      Z_dense(i, j) = Z_dense(i, j)+1/sqrt(eval(k))*S_dense(i, k)*S_dense(j, k)
    end forall
  enddo
  !call print_matrix_python(X_dense, "X (S^{-1/2})")
  !call print_matrix_python(Y_dense, "Y (S^{1/2})")
  !call print_matrix_python(matmul(Y_dense, Y_dense), "Y2 (S)")
  !call print_matrix_python(matmul(X_dense, Y_dense), "X*Y (I)")

  Y_norm = 0
  max_diff = 0
  do i = 1, N
    do j = 1, N
      abs_diff = abs(Y_dense(i, j)-get(Y, i, j))
      Y_norm = Y_norm+abs_diff**2
      if(abs_diff > max_diff) then
        max_diff = abs_diff
      endif
    end do
  end do
  Y_norm = Y_norm/(N*N)
  write(*, "(A,ES10.3)") "Y max diff = ", max_diff
  write(*, "(A,ES10.3)") "||Y-Y_dense||_{F}/N**2 = ", Y_norm

  Z_norm = 0
  max_diff = 0
  do i = 1, N
    do j = 1, N
      abs_diff = abs(Y_dense(i, j)-get(Y, i, j))
      Z_norm = Z_norm+abs_diff**2
      if(abs_diff > max_diff) then
        max_diff = abs_diff
      endif
    end do
  end do
  Z_norm = Z_norm/(N*N)
  write(*, "(A,ES10.3)") "Z max diff = ", max_diff
  write(*, "(A,ES10.3)") "||Y-Y_dense||_{F}/N**2 = ", Z_norm

  if(Y_norm > 1e-8_spamm_kind .or. Z_norm > 1e-8_spamm_kind) then
    write(*, "(A)") "mismatch"
    error stop
  else
    write(*, "(A)") "test passed"
  endif
#endif

  call delete(S)
  call delete(Y)
  call delete(Z)

end program test
