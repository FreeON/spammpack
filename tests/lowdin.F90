program test

  use spammpack
  use test_utilities

  implicit none

  interface
    subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LIWORK, LWORK, N
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
    end subroutine dsyevd
  end interface

  integer, parameter :: N = 10
  integer, parameter :: LWORK = 1+6*N+2*N**2
  integer, parameter :: LIWORK = 3+5*N

  type(spamm_matrix_2nd_order), pointer :: S, X
  real(kind(0d0)), dimension(N, N) :: S_dense, X_dense, Y_dense
  real(kind(0d0)), dimension(N) :: eval
  real(kind(0d0)), dimension(N, N) :: evec
  real(kind(0d0)), dimension(LWORK) :: work
  integer, dimension(LIWORK) :: iwork
  integer :: info
  integer :: i, j, k
  real(kind(0d0)) :: abs_diff, max_diff

  call random_number(eval)
  call random_number(evec)

  eval = eval/2+0.4
  write(*, "(A,ES10.3)") "condition number = ", maxval(eval)/minval(eval)

  do k = 1, N
    evec(k, :) = evec(k, :)/sqrt(sum(evec(k, :)*evec(k,:)))
  enddo

  call print_matrix_python(reshape(eval, [ 1, size(eval) ]), "lambda")
  call print_matrix_python(evec, "basis")

  S_dense = 0
  do k = 1, N
    forall(i = 1:N, j = 1:N)
      S_dense(i, j) = S_dense(i, j)+eval(k)*evec(k, i)*evec(k, j)
    end forall
  enddo

  call print_matrix_python(S_dense, "S")

  S => spamm_convert_dense_to_matrix_2nd_order(S_dense)
  X => spamm_inverse_sqrt_schulz(S)

#ifdef LAPACK_FOUND
  call dsyevd("V", "U", N, S_dense, N, eval, work, LWORK, iwork, LIWORK, info)

  write(*, "(A)") "condition number: "//to_string(eval(N)/eval(1))

  call print_matrix_python(reshape(eval, [ 1, size(eval) ]), "eigenvalue")
  call print_matrix_python(S_dense, "eigenvector")

  X_dense = 0
  Y_dense = 0
  do k = 1, N
    forall(i = 1:N, j = 1:N)
      X_dense(i, j) = X_dense(i, j)+1/sqrt(eval(k))*S_dense(i, k)*S_dense(j, k)
      Y_dense(i, j) = Y_dense(i, j)+sqrt(eval(k))*S_dense(i, k)*S_dense(j, k)
    end forall
  enddo
  call print_matrix_python(X_dense, "X (S^{-1/2})")
  call print_matrix_python(Y_dense, "Y (S^{1/2})")
  call print_matrix_python(matmul(Y_dense, Y_dense), "Y2 (S)")
  call print_matrix_python(matmul(X_dense, Y_dense), "X*Y (I)")

  max_diff = 0
  do i = 1, N
    do j = 1, N
      abs_diff = abs(X_dense(i, j)-get(X, i, j))
      if(abs_diff > max_diff) then
        max_diff = abs_diff
      endif
    end do
  end do
  write(*, "(A,ES10.3)") "max diff = ", max_diff
  if(max_diff > 1d-10) then
    write(*, "(A)") "mismatch"
    error stop
  endif
#endif

  call delete(S)
  call delete(X)

end program test
