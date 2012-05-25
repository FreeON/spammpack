program spamm_multiply

  use spammpack
  use spammtests

  implicit none

  integer :: N
  integer :: testresult = 0
  integer, parameter :: NThreads = 4

  real(SpAMM_DOUBLE) :: max_norm

  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: A_dense
  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: B_dense
  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: C_dense

  real(SpAMM_DOUBLE), dimension(:,:), allocatable :: C_diff

  real(SpAMM_KIND), dimension(:,:), allocatable :: A_dense_padded
  real(SpAMM_KIND), dimension(:,:), allocatable :: B_dense_padded
  real(SpAMM_KIND), dimension(:,:), allocatable :: C_dense_padded

  type(QuTree), pointer :: A => null()
  type(QuTree), pointer :: B => null()
  type(QuTree), pointer :: C => null()
  type(QuTree), pointer :: C_reference => null()

  integer :: i, j

  call load_matrix("testmatrix_2.coor", A_dense)
  call load_matrix("testmatrix_2.coor", B_dense)

  N = size(A_dense, 1)
  allocate(C_dense(N, N))
  C_dense = SpAMM_ZERO

  call SpAMM_Init_Globals(N, NThreads)

  allocate(A_dense_padded(N, N))
  allocate(B_dense_padded(N, N))
  allocate(C_dense_padded(N, N))

  A_dense_padded = SpAMM_ZERO
  B_dense_padded = SpAMM_ZERO
  C_dense_padded = SpAMM_ZERO

  A_dense_padded = A_dense
  B_dense_padded = B_dense

  A => SpAMM_Convert_Dense_2_QuTree(A_dense_padded)
  B => SpAMM_Convert_Dense_2_QuTree(B_dense_padded)
  call New(C)

  !$OMP PARALLEL
  !$OMP SINGLE

  call Multiply(A, B, C)

  !$OMP END SINGLE
  !$OMP END PARALLEL

  C%Norm = SQRT(Norm(C))

  CALL SpAMM_Time_Stamp()

  C_dense = matmul(A_dense, B_dense)
  C_dense_padded = C_dense
  C_reference => SpAMM_Convert_Dense_2_QuTree(C_dense_padded)

  call Add(C, C_reference, -SpAMM_ONE, SpAMM_ONE)

  max_norm = SQRT(Norm(C))

  write(*, *) "diff norm (C) = ", max_norm

  ! Exit with some error code.
  call spamm_exit(testresult)

end program spamm_multiply
