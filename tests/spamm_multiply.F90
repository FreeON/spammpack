#ifndef TEST_REPEAT
#define TEST_REPEAT 1
#endif

program spamm_multiply

  use spammpack
  use spammtests

  implicit none

  integer :: N
  integer :: test_repeat
  integer :: testresult = 0

#ifdef _OPENMP
  integer :: min_threads
  integer :: max_threads
  integer :: num_threads
#endif

  type(SpAMM_Norm) :: norms

  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: A_dense
  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: B_dense
  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: C_dense

  real(SpAMM_KIND), dimension(:,:), allocatable :: A_dense_padded
  real(SpAMM_KIND), dimension(:,:), allocatable :: B_dense_padded
  real(SpAMM_KIND), dimension(:,:), allocatable :: C_dense_padded

  type(QuTree), pointer :: A => null()
  type(QuTree), pointer :: B => null()
  type(QuTree), pointer :: C => null()
  type(QuTree), pointer :: C_reference => null()

  character(len = 1000) :: inputbuffer
  character(len = 1000) :: matrixfilename

  call get_command_argument(1, matrixfilename)

  if(matrixfilename == "") then
    matrixfilename = "testmatrix_random_1024x1024.coor"
  endif

#ifdef _OPENMP
  call get_command_argument(2, inputbuffer)
  read(inputbuffer, "(I)") num_threads
#endif

  call load_matrix(matrixfilename, A_dense)
  call load_matrix(matrixfilename, B_dense)

  N = size(A_dense, 1)
  allocate(C_dense(N, N))
  C_dense = SpAMM_ZERO

  write(*, *) "read matrix N = ", N

  ! Get new, padded matrix size.
#ifdef _OPENMP
  call SpAMM_Init_Globals(N, num_threads)
#else
  call SpAMM_Init_Globals(N)
#endif

  write(*, *) "padded matrix to N = ", N

  allocate(A_dense_padded(N, N))
  allocate(B_dense_padded(N, N))
  allocate(C_dense_padded(N, N))

  A_dense_padded = SpAMM_ZERO
  B_dense_padded = SpAMM_ZERO
  C_dense_padded = SpAMM_ZERO

  A_dense_padded = A_dense
  B_dense_padded = B_dense

  write(*, *) "converting matrices to quadtree"
  A => SpAMM_Convert_Dense_2_QuTree(A_dense_padded)
  B => SpAMM_Convert_Dense_2_QuTree(B_dense_padded)
  call New(C)

#if defined(_OPENMP)
  if(num_threads == 0) then
    min_threads = 1
    max_threads = omp_get_num_procs()
  else
    min_threads = num_threads
    max_threads = num_threads
  endif

  write(*, "(A,I2,A,I2)") "cycling number of threads between ", min_threads, " and ", max_threads

  do num_threads = min_threads, max_threads
    CALL SpAMM_Set_Num_Threads(num_threads)
    CALL SpAMM_Timer_Reset()
#endif

    write(*, "(A,I4)") "repeat multiply ", TEST_REPEAT
    do test_repeat = 1, TEST_REPEAT
      call Multiply(A, B, C, LocalThreshold = 1e-5)
    enddo

    CALL SpAMM_Time_Stamp()

#ifdef VERIFY_RESULT
    C_dense = matmul(A_dense, B_dense)
    C_dense_padded = C_dense
    C_reference => SpAMM_Convert_Dense_2_QuTree(C_dense_padded)

    call Add(C, C_reference, -SpAMM_ONE, SpAMM_ONE)

    norms = Norm(A)
    write(*, "(A,F22.12)") "F-norm (A)   = ", sqrt(norms%FrobeniusNorm)
    write(*, "(A,F22.12)") "max-norm (A) = ", norms%MaxNorm

    norms = Norm(C)

    write(*, "(A,F22.12)") "F-norm (C)   = ", sqrt(norms%FrobeniusNorm)
    write(*, "(A,F22.12)") "max-norm (C) = ", norms%MaxNorm
#endif

#if defined(_OPENMP)
  enddo
#endif

  ! Exit with some error code.
  call spamm_exit(testresult)

end program spamm_multiply
