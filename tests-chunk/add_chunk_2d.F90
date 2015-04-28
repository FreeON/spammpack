program test

  use spamm_chunk

  implicit none

  type(chunk_2d_t), pointer :: A, B, C
  double precision :: A_dense(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)
  double precision :: B_dense(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)
  double precision :: C_dense(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)
  integer :: i, j

  allocate(A)
  allocate(B)
  call init_2d_random(A)
  call init_2d_random(B)

  nullify(C)
  C => chunk_add_2d_2d(A, B)

  A_dense = chunk_2d_to_dense(A)
  B_dense = chunk_2d_to_dense(B)

  do i = 1, SPAMM_CHUNK_SIZE
     do j = 1, SPAMM_CHUNK_SIZE
        if(abs(A_dense(i, j)-B_dense(i, j)) > 1e-12) then
           write(*, *) "matrix element mismatch"
           write(*, *) trim(chunk_2d_to_string(A))
           write(*, *) trim(chunk_2d_to_string(B))
           write(*, *) A_dense
           write(*, *) B_dense
           error stop
        endif
     enddo
  enddo

  write(*, *) "matrices match"

end program test
