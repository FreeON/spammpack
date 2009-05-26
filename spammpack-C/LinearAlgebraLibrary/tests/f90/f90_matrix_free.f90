program matrix_free

  use lal_types

  type(lal_matrix) :: A

  if(f90_lal_allocate(10, 10, A) /= 0) then
    call exit(-1)
  endif

  call f90_lal_zero(A)

  call f90_lal_free(A)

end program matrix_free
