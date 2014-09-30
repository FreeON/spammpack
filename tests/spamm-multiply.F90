#ifndef TEST_REPEAT
#define TEST_REPEAT 1
#endif

#undef VERIFY_RESULT

program spamm_multiply

  use spammpack
  use test_utilities

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  real(spamm_kind), parameter :: tolerance = 1d-8

  integer :: M, N
  integer :: test_repeat
  integer :: testresult = 0

  real(kind(0d0)), dimension(:,:), allocatable :: A_dense

  type(spamm_matrix_2nd_order), pointer :: A => null()
  type(spamm_matrix_2nd_order), pointer :: B => null()
  type(spamm_matrix_2nd_order), pointer :: C => null()

#ifdef VERIFY_RESULT
  real(spamm_kind) :: norms
  real(kind(0d0)), dimension(:,:), allocatable :: C_dense
  type(spamm_matrix_2nd_order), pointer :: C_reference => null()
#endif

#ifdef _OPENMP
  integer :: min_threads
  integer :: max_threads
  integer :: num_threads = 1
  character(len = 1000) :: inputbuffer
#endif
  character(len = 1000) :: matrixfilename

  if(command_argument_count() > 0) then
    call get_command_argument(1, matrixfilename)
  else
    write(*, *) "Please specify a matrix file (in MM format)"
    error stop
  endif

#ifdef _OPENMP
  if(command_argument_count() > 1) then
    call get_command_argument(2, inputbuffer)
    read(inputbuffer, "(I3)") num_threads
  else
    num_threads = omp_get_max_threads()
  endif

  write(*, *) "setting num_threads to ", num_threads
#endif

  call read_MM(matrixfilename, A_dense)
  !call load_matrix(matrixfilename, A_dense)
  !call load_matrix_binary(matrixfilename, B_dense)

  M = size(A_dense, 1)
  N = size(A_dense, 2)

  write(*, *) "read matrix N = ", N, ", M = ", M

  write(*, *) "converting matrices to quadtree"
  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => A
  write(*, *) "done converting"

#if defined(_OPENMP)
  if(num_threads == 0) then
    min_threads = 1
    max_threads = omp_get_num_procs()
  else
    min_threads = num_threads
    max_threads = num_threads
  endif

  write(*, "(A,I3,A,I3)") "cycling number of threads between ", min_threads, " and ", max_threads

  do num_threads = min_threads, max_threads
    CALL omp_set_num_threads(num_threads)
    CALL SpAMM_Timer_Reset()
#endif

    write(*, "(A,I4)") "repeat multiply ", TEST_REPEAT
    do test_repeat = 1, TEST_REPEAT
      call multiply(A, B, C, tolerance)
      write(*, "(A,F14.1)") "operation count = ", C%number_operations
      write(*, "(A,F14.4)") "count/N         = ", C%number_operations/dble(N)
      write(*, "(A,F14.4)") "count/N^3       = ", C%number_operations/dble(N)**3
    enddo

    !CALL SpAMM_Time_Stamp()

#ifdef VERIFY_RESULT
    allocate(C_dense(M, N))
    C_dense = matmul(A_dense, A_dense)
    C_reference => spamm_convert_dense_to_matrix_2nd_order(C_dense)

    norms = Norm(A)
    write(*, "(A,F22.12)") "F-norm (A)             = ", sqrt(norms)

    norms = Norm(B)
    write(*, "(A,F22.12)") "F-norm (B)             = ", sqrt(norms)

    norms = Norm(C)
    write(*, "(A,F22.12)") "F-norm (C)             = ", sqrt(norms)

    norms = Norm(C_reference)
    write(*, "(A,F22.12)") "F-norm (C_reference)   = ", sqrt(norms)

    call Add(C, C_reference, -SpAMM_ONE, SpAMM_ONE)

    norms = Norm(C)
    write(*, "(A,F22.12)") "F-norm (diff)          = ", sqrt(norms)
#endif

#if defined(_OPENMP)
  enddo
#endif

  ! Exit with some error code.
  call spamm_exit(testresult)

end program spamm_multiply
