module lal

  type lal_matrix
    character(len=4) :: address
  end type lal_matrix

  contains

  function lal_get (i, j, A)
    integer          :: i, j
    type(lal_matrix) :: A
  end function lal_get

  subroutine lal_allocate (i, j, A)
    integer          :: i, j
    type(lal_matrix) :: A
  end subroutine lal_allocate

end module lal
