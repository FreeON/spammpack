program spamm_multiply

  use spammpack
  use test_utilities

  implicit none

  integer :: N
  integer, parameter :: NThreads = 8

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

  !write(*, *) "N = ", N

  call SpAMM_Init_Globals(N, NThreads)
  !write(*, *) "N = ", N

  allocate(A_dense_padded(N, N))
  allocate(B_dense_padded(N, N))
  allocate(C_dense_padded(N, N))

  A_dense_padded = SpAMM_ZERO
  B_dense_padded = SpAMM_ZERO
  C_dense_padded = SpAMM_ZERO

  A_dense_padded = A_dense
  B_dense_padded = B_dense

  !write(*,*) A_dense(1, 1), B_dense(1, 1)
  !write(*,*) A_dense_padded(1, 1), B_dense_padded(1, 1)

  !call dgemm('N', 'N', N, N, N, 1.0d0, A_dense, N, B_dense, N, 1.0d0, C_dense, N)
  !C_dense = matmul(A_dense, B_dense)
  !call sgemm('N', 'N', N, N, N, 1.0, A_dense_padded, N, B_dense_padded, N, 1.0, C_dense_padded, N)
  !C_dense_padded = matmul(A_dense_padded, B_dense_padded)

  !allocate(C_diff(N, N))
  !write(*,*) C_dense(1, 1), C_dense_padded(1, 1)
  !C_diff = C_dense-C_dense_padded

  !max_norm = 0.0
  !do i = 1, N
  !  do j = 1, N
  !    max_norm = max(max_norm, C_diff(i, j))
  !  enddo
  !enddo

  !write(*,*) "max_norm = ", max_norm
  !stop

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

  !write(*, *) "here"
  C_dense = matmul(A_dense, B_dense)
  !write(*, *) "here 2"
  C_dense_padded = C_dense
  !write(*, *) "here 3"
  C_reference => SpAMM_Convert_Dense_2_QuTree(C_dense_padded)
  !write(*, *) "here 4"

  call Add(C, C_reference, -SpAMM_ONE, SpAMM_ONE)

  max_norm = SQRT(Norm(C))

  write(*, *) "diff norm (C) = ", max_norm

end program spamm_multiply
