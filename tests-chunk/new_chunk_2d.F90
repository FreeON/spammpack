program test

  use spammchunk

  implicit none

  type(chunk_2d), pointer :: A

  allocate(A)

  write(*, "(A)") trim(chunk_2d_to_string(A))

end program test
