program test

  character(len = 10) :: a
  character(len = :), allocatable :: b

  a = "12345"
  b = trim(a)
  write(*, "(A)") "'"//a//"' '"//b//"'"
  write(*, *) len(a), len(b)

  if(len(b) /= len_trim(a)) then
    write(*, *) "wrong string length"
    error stop
  endif

end program test
