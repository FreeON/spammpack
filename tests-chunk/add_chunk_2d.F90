program test

  use spammchunk

  implicit none

  type(chunk_2d), pointer :: A, B
  real(kind(0d0)) :: A_dense(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)
  real(kind(0d0)) :: B_dense(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)

  allocate(A)
  call init_2d_random(A)

  B => chunk_add_2d_2d(A, A)

end program test
