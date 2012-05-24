program multiply

  use spammpack
  use load_matrix

  implicit none

  integer :: N
  real(kind = 8), dimension(:,:), allocatable :: A
  real(kind = 8), dimension(:,:), allocatable :: B
  real(kind = 8), dimension(:,:), allocatable :: C

  call load("testmatrix_1.coor", A)
  call load("testmatrix_1.coor", B)
  call load("testmatrix_1.coor", C)

  N = size(A, 1)

  call SpAMM_Init_Globals(N)

end program multiply
