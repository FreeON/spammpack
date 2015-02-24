program test

  use spammchunk

  type(chunk_2d), pointer :: A
  real(kind(0d0)) :: B(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)

  allocate(A)

  call init_2d_random(A)
  B = chunk_2d_to_dense(A)

end program test
