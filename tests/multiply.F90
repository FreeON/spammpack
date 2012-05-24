program multiply

  use spammpack
  use test_utilities

  implicit none

  integer :: N

  real(kind = 8) :: alpha
  real(kind = 8) :: beta

  real(kind = 8), dimension(:,:), allocatable :: A
  real(kind = 8), dimension(:,:), allocatable :: B
  real(kind = 8), dimension(:,:), allocatable :: C

  alpha = 1.2
  beta = 0.5

  call load_matrix("testmatrix_1.coor", A)
  call load_matrix("testmatrix_1.coor", B)
  call load_matrix("testmatrix_1.coor", C)

  N = size(A, 1)

  call SpAMM_Init_Globals(N)

end program multiply
