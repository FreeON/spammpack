program test

  use spammpack

  implicit none

  type(chunk_2d), pointer :: A

  A => new_chunk_2d(128)
  write(*, "(A)") trim(chunk_2d_to_string(A))
  write(*, "(A)") trim(to_string(storage_size(A)/8))//" bytes"

end program test
