program test

  use spammchunk

  type(chunk_2d), pointer :: A
  real(kind(0d0)) :: B(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)

  allocate(A)

  call init_2d_random(A)
  write(*, "(A)") trim(chunk_2d_to_string(A))

  B = chunk_2d_to_dense(A)

  do i = 1, SPAMM_CHUNK_SIZE
     do j = 1, SPAMM_CHUNK_SIZE
        if(B(i, j) /= chunk_2d_get(i, j, A)) then
           write(*, *) "matrix element mismatch at ", i, j
           write(*, *) "chunk = ", chunk_2d_get(i, j, A)
           write(*, *) "dense = ", B(i, j)
           error stop
        endif
     enddo
  enddo

end program test
