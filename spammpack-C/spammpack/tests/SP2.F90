function IdmpCnvrgChck (Occ0, Occ1, Occ2, Occ3, iMin, i) result(converged)

  use spammpack

  logical                        :: converged
  real(SpAMM_DOUBLE), intent(in) :: Occ0, Occ1, Occ2, Occ3
  integer, intent(in)            :: iMin, i

  real(SpAMM_DOUBLE) :: IdmpErrCurr, IdmpErrOld

  converged = .false.
  if(i >= iMin) then
    IdmpErrCurr = abs(Occ0-Occ1) ! |Tr(P*Q)|
    if(IdmpErrCurr < 0.01D0) then
      IdmpErrOld = abs(Occ2-Occ3)
      if(IdmpErrCurr >= IdmpErrOld) then
        converged = .true.
      endif
    endif
  endif

end function IdmpCnvrgChck

program spamm_SP2

  use spammpack
  use spammtests

  implicit none

  logical :: IdmpCnvrgChck

  integer :: N, N_padded
  integer :: i, j
  integer :: test_repeat
  integer :: testresult = 0

#ifdef _OPENMP
  integer :: min_threads
  integer :: max_threads
  integer :: num_threads
#endif

  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: P_dense

  real(SpAMM_KIND), dimension(:,:), allocatable :: P_dense_padded

  type(SpAMM_Matrix) :: P
  type(SpAMM_Matrix) :: P2

  character(len = 1000) :: inputbuffer
  character(len = 1000) :: matrixfilename

  integer :: intbuffer

  integer :: I2, MM, iMin
  real(SpAMM_KIND) :: TrE, HalfNe, Nocc
  real(SpAMM_DOUBLE) :: Ne, Lambda, idempotency_error, Occ0, Occ1, Occ2, Occ3

  call get_command_argument(1, matrixfilename)

  if(matrixfilename == "") then
    matrixfilename = "testmatrix_random_1024x1024.coor"
  endif

  call get_command_argument(2, inputbuffer)
  read(inputbuffer, "(I8)"), intbuffer
  Ne = intbuffer

  if(Ne <= 0) then
    write(*, *) "Number of electrons missing or wrong"
    call SpAMM_Exit(1)
  else
    write(*, "(A,D12.2)") "Setting Ne = ", Ne
  endif

#ifdef _OPENMP
  call get_command_argument(3, inputbuffer)
  read(inputbuffer, "(I2)") num_threads
#endif

  call load_matrix(matrixfilename, P_dense)

  N = size(P_dense, 1)

  write(*, *) "read matrix N = ", N

  ! Get new, padded matrix size.
  N_padded = N
#ifdef _OPENMP
  call SpAMM_Init_Globals(N_padded, num_threads)
#else
  call SpAMM_Init_Globals(N_padded)
#endif

  write(*, *) "padded matrix to N_padded = ", N_padded

  allocate(P_dense_padded(N_padded, N_padded))

  P_dense_padded = SpAMM_ZERO
  P_dense_padded = P_dense

  write(*, *) "converting matrices to quadtree"
  P = SpAMM_Convert_Dense_to_SpAMM(P_dense_padded)

  ! Allocate new P2.
  call New(P2)

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
    call SpAMM_Set_Num_Threads(num_threads)
    call SpAMM_Timer_Reset()
#endif

    HalfNe = SpAMM_Half*Ne
    Occ0 = SpaMM_Zero
    Occ1 = SpaMM_Zero
    Occ2 = SpaMM_Zero
    Occ3 = SpaMM_Zero
    iMin = 20

    write(*, *) "          i   occupation      occError"
    DO i = 1, 100
      call SpAMM_TC2(P, P2, HalfNe, Nocc)
      write(*, *) i, ABS(Nocc*SpAMM_Two), ABS(Nocc*SpAMM_Two-Ne)
      Occ0=Nocc
      IF(IdmpCnvrgChck(Occ0, Occ1, Occ2, Occ3, iMin, i)) THEN
        write(*, *) "converged in ", i, "iterations"
        write(*, *) "Idempotency error = ", abs(Occ0-Occ1)
        write(*, *) "Previous idempotency error = ", abs(Occ2-Occ3)
        exit
      endif
      Occ3 = Occ2
      Occ2 = Occ1
      Occ1 = Occ0
    enddo

    CALL SpAMM_Time_Stamp()

#if defined(_OPENMP)
  enddo
#endif

  ! Exit with some error code.
  call spamm_exit(testresult)

end program spamm_SP2
