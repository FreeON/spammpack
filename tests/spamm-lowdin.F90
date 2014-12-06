program spamm_lowdin

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

  real(kind(0d0)), parameter :: tolerance = 1d-8

  type(spamm_matrix_2nd_order), pointer :: S => null(), X => null(), &
    Y => null(), Z => null(), Id => null()
  real(kind(0d0)), allocatable :: S_dense(:, :), Id_dense(:, :)
  character(len = 1000) :: matrix_filename
  real :: start_time, end_time

  real(kind(0d0)) :: max_diff

#ifdef LAPACK_FOUND
  integer :: N
  integer :: LWORK
  integer :: LIWORK

  real(kind(0d0)), allocatable :: eval(:), evec(:, :), work(:)
  integer, allocatable :: iwork(:)
  integer :: info
#endif

  if(command_argument_count() == 0) then
    write(*, "(A)") "Please specify the matrix file to use (the overlap matrix)"
    error stop
  else if(command_argument_count() > 1) then
    write(*, "(A)") "Too many command line arguments. One is plenty."
    error stop
  else
    call get_command_argument(1, matrix_filename)
  endif

  call read_MM(matrix_filename, S_dense)
  S => spamm_convert_dense_to_matrix_2nd_order(S_dense)

  call cpu_time(start_time)
  call spamm_inverse_sqrt_schulz(S, Y, Z, tolerance)
  call cpu_time(end_time)

  write(*, "(A)") "Y (S^{+1/2}) fillin: "//to_string(Y%number_nonzeros)
  write(*, "(A)") "Z (S^{-1/2}) fillin: "//to_string(Z%number_nonzeros)
  write(*, "(A)") "CPU time: "//to_string(end_time-start_time)

  Id => spamm_identity_matrix(S%M, S%N)
  call multiply(Y, Z, X)
  call add(X, Id, +1.0d0, -1.0d0)
  max_diff = absmax(X)
  write(*, "(A)") "|Id-S^{-1/2} S^{1/2}|_{max}: "//to_string(max_diff)

#ifdef LAPACK_FOUND
  N = size(S_dense, 1)
  LWORK = 1+6*N+2*N**2
  LIWORK = 3+5*N

  allocate(eval(N))
  allocate(evec(N, N))
  allocate(work(LWORK))
  allocate(iwork(LIWORK))

  write(*, "(A)") "Getting condition number..."
  call dsyevd("V", "U", N, S_dense, N, eval, work, LWORK, iwork, LIWORK, info)

  write(*, "(A)") "Condition number: "//to_string(eval(N)/eval(1))

  deallocate(iwork)
  deallocate(work)
  deallocate(evec)
  deallocate(eval)
#endif

end program spamm_lowdin
