program set_element

  use lal

  implicit none

  type(lal_matrix) :: A

  if (f90_lal_allocate(10, 10, A) /= 0) then
    call exit(-1)
  endif

  call f90_lal_zero(A);
  call f90_lal_set(1, 2, 5d0, A);

  if(f90_lal_get(1, 2, A) /= 5d0) THEN
    call exit(-1)
  endif

end program set_element
