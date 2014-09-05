program test

  ! This is a long line (longer than 132 characters).
  character(len = *), parameter :: long_line_string = "This is the really long line. The length is longer than 132 characters to test how forgiving the compiler is with this kind of line length."

  write(*, *) long_line_string

end program test
