program matrix_allocation

  use lal

  implicit none

  integer          :: i, j
  type(lal_matrix) :: A

  i = 10
  j = 10

  if(f90_lal_allocate(i, j, A) /= 0) then
    call exit(-1)
  endif

end program
