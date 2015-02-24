program test

  use spammchunk

  implicit none

  type(chunk_2d), pointer :: A

  allocate(A)

  write(*, "(A)") trim(chunk_2d_to_string(A))
  write(*, "(A)") "total: "//trim(to_string(storage_size(A)/8))//" bytes"

end program test
