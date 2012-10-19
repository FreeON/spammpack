program fortran_C

  use iso_c_binding

  integer(kind = c_int) :: n
  integer(kind = c_int64_t) :: i

  n = c_func_1()
  write(*,*) "n = ", n

  call c_func_2(n)
  write(*,*) "n = ", n

  i = c_func_3()
  write(*,*) "i = ", i

  call c_func_4(n, i)
  write(*,*) "i = ", i

  call c_func_5(i)
  write(*,*) "n = ", n

end program fortran_C
