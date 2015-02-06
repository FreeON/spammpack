program spamm_lowdin

  use spammpack
  use test_utilities

  implicit none

#ifdef LAPACK_FOUND
  interface
     subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
       character          JOBZ, UPLO
       integer            INFO, LDA, LIWORK, LWORK, N
       integer            IWORK( * )
       double precision   A( LDA, * ), W( * ), WORK( * )
     end subroutine dsyevd
  end interface
#endif

  ! COND(5)
  real(kind(0d0)), parameter :: tolerance = 1d-8

  type(tree_2d_symmetric), pointer :: S => null()
  type(tree_2d_symmetric), pointer :: X => null()
  type(tree_2d_symmetric), pointer :: Y => null()
  type(tree_2d_symmetric), pointer :: Z => null()
  type(tree_2d_symmetric), pointer :: Z2 => null()
  type(tree_2d_symmetric), pointer :: Id => null()

  real(kind(0d0)), allocatable :: S_dense(:, :)
  character(len = 1000) :: matrix_filename
  real :: start_time, end_time
  integer :: i
  real(kind(0d0)) :: max_diff
  logical :: input_file_exists

#ifdef LAPACK_FOUND
  integer :: N
  integer :: LWORK
  integer :: LIWORK

  real(kind(0d0)), allocatable :: eval(:)
  real(kind(0d0)), allocatable :: SHalf(:,:)
  real(kind(0d0)), allocatable :: SHlfI(:,:)
  real(kind(0d0)), allocatable :: evec(:, :)
  real(kind(0d0)), allocatable :: work(:)

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

  inquire(file=matrix_filename, exist=input_file_exists)
  if(.not. input_file_exists) then
     write(*, "(A)") "Input matrix file does not exist"
     error stop
  endif

  call read_MM(matrix_filename, S_dense)
  S => spamm_convert_dense_to_matrix_2nd_order(S_dense)

  write(*,*)' computing inverse now .... '
  write(*,*)' using spamm threshold = ',to_string(tolerance)

  call cpu_time(start_time)
  call spamm_inverse_sqrt_schulz(S, Y, Z, tolerance)  !, schulz_threshold=1.d-12)
  call cpu_time(end_time)

  write(*, "(A)") "Y (S^{+1/2}) fillin: "//to_string(Y%decoration%number_nonzeros)
  write(*, "(A)") "Z (S^{-1/2}) fillin: "//to_string(Z%decoration%number_nonzeros)
  write(*, "(A)") "Schulz CPU time: "//to_string(end_time-start_time)

#ifdef LAPACK_FOUND
  N = size(S_dense, 1)
  LWORK = 1+6*N+2*N**2
  LIWORK = 3+5*N

  allocate(eval(N))
  allocate(work(LWORK))
  allocate(iwork(LIWORK))

  write(*, "(A)") "Getting condition number..."

  call cpu_time(start_time)
  call dsyevd("V", "U", N, S_dense, N, eval, work, LWORK, iwork, LIWORK, info)
  call cpu_time(end_time)

  deallocate(iwork)
  deallocate(work)
  allocate(evec(N, N))
  allocate(SHalf(N, N))
  allocate(SHlfI(N, N))

  write(*,*)' MINMAX =',Eval(1),Eval(N)

  evec = S_dense

  do i = 1, N
     S_dense(:,I) = S_dense(:,i)*sqrt(eval(i))
  enddo

  SHalf = matmul(evec,transpose(S_dense))

  do i = 1, N
     S_dense(:,I) = S_dense(:,i)/eval(i)
  enddo

  SHlfI = matmul(evec,transpose(S_dense))

  deallocate(evec)
  deallocate(eval)
  deallocate(S_dense)

  Y  => spamm_convert_dense_to_matrix_2nd_order(SHalf)
  Z2 => spamm_convert_dense_to_matrix_2nd_order(SHlfI)
#endif

  write(*,*)' checking with Z[SpAMM] '
  Id => spamm_identity_matrix(S%decoration%M, S%decoration%N)
  call multiply(Y, Z, X , 0D0)
  call add(X, Id, +1.0d0, -1.0d0)
  max_diff = absmax(X)
  write(*, "(A)") "|Id-S^{-1/2} S^{1/2}_dsyev|_{max}: "//to_string(max_diff)

  write(*,*)' checking with Z[DSYEV] '
  Id => spamm_identity_matrix(S%decoration%M, S%decoration%N)
  call multiply(Y, Z2, X , 0D0)
  call add(X, Id, +1.0d0, -1.0d0)
  max_diff = absmax(X)
  write(*, "(A)") "|Id-S^{-1/2}_dsyev S^{1/2}_dsyev|_{max}: "//to_string(max_diff)

end program spamm_lowdin
