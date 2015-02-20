program test

  use spammpack

  implicit none

  type(chunk_2d), pointer :: A

  A => new_chunk_2d(128)
  write(*, "(A)") trim(A%to_string())
  write(*, "(A)") trim(to_string(storage_size(A)/8))//" bytes"

end program test
