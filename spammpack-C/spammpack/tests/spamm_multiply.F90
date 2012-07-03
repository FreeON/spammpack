#ifndef TEST_REPEAT
#define TEST_REPEAT 1
#endif

#define VERIFY_RESULT

program spamm_multiply

  use spammpack
  use spammtests

  implicit none

  integer :: N, N_padded
  integer :: test_repeat
  integer :: testresult = 0

  real(SpAMM_SINGLE), dimension(:,:), allocatable :: A_dense
  real(SpAMM_SINGLE), dimension(:,:), allocatable :: B_dense
  real(SpAMM_SINGLE), dimension(:,:), allocatable :: C_dense

  type(SpAMM_Matrix) :: A
  type(SpAMM_Matrix) :: B
  type(SpAMM_Matrix) :: C

#ifdef VERIFY_RESULT
  !type(SpAMM_Norm) :: norms
  integer :: i, j
  type(SpAMM_Matrix) :: C_reference
#endif

#ifdef _OPENMP
  integer :: min_threads
  integer :: max_threads
  integer :: num_threads
  character(len = 1000) :: inputbuffer
#endif
  character(len = 1000) :: matrixfilename

  call get_command_argument(1, matrixfilename)

  if(matrixfilename == "") then
    matrixfilename = "testmatrix_random_1024x1024.coor"
  endif

#ifdef _OPENMP
  call get_command_argument(2, inputbuffer)
  read(inputbuffer, "(I3)") num_threads
#endif

  call load_matrix(matrixfilename, A_dense)
  !call load_matrix_binary(matrixfilename, A_dense)
  !call load_matrix_binary(matrixfilename, B_dense)

  N = size(A_dense, 1)
  allocate(B_dense(N, N))
  B_dense = A_dense
  allocate(C_dense(N, N))
  C_dense = 0.0D0

  write(*, *) "read matrix N = ", N

  write(*, *) "converting matrices to quadtree"
  A = SpAMM_Convert_Dense_to_SpAMM(A_dense)
  B = SpAMM_Convert_Dense_to_SpAMM(B_dense)
  call New(size(A_dense, 1), size(B_dense, 2), C)
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
    CALL SpAMM_Set_Num_Threads(num_threads)
    !CALL SpAMM_Timer_Reset()
#endif

    write(*, "(A,I4)") "repeat multiply ", TEST_REPEAT
    do test_repeat = 1, TEST_REPEAT
      call Multiply(A, B, C, tolerance = 1e-7)
    enddo

    !CALL SpAMM_Time_Stamp()

#ifdef VERIFY_RESULT
    C_dense = matmul(A_dense, B_dense)
    C_reference = SpAMM_Convert_Dense_to_SpAMM(C_dense)

    write(*, "(A,F22.12)") "F-norm (A)             = ", SpAMM_Norm(A)
    write(*, "(A,F22.12)") "F-norm (B)             = ", SpAMM_Norm(B)
    write(*, "(A,F22.12)") "F-norm (C)             = ", SpAMM_Norm(C)
    write(*, "(A,F22.12)") "F-norm (C_reference)   = ", SpAMM_Norm(C_reference)

    call Add(C, C_reference, -SpAMM_ONE, SpAMM_ONE)

    write(*, "(A,F22.12)") "F-norm (diff)          = ", SpAMM_Norm(C)
#endif

#if defined(_OPENMP)
  enddo
#endif

  ! Exit with some error code.
  call spamm_exit(testresult)

end program spamm_multiply
