program test

  use spamm_chunk

  implicit none

  type(chunk_2d_t), pointer :: A

  allocate(A)

  write(*, "(A)") trim(chunk_2d_to_string(A))
  write(*, "(A)") trim(chunk_2d_memory_layout(A))

end program test
