program matrix_allocation

  use lal

  implicit none

  integer          :: i, j
  type(lal_matrix) :: A

  i = 10
  j = 10

  call lal_allocate(i, j, A)

end program
