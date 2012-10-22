program chunk_test

  use :: spammpack
  use, intrinsic :: iso_c_binding

  integer :: number_dimensions = 2
  integer :: N_contiguous = 8
  integer, dimension(:), pointer :: N_lower
  integer, dimension(:), pointer :: N_upper
  real*4, dimension(:,:), pointer :: A, B, C

  type(c_ptr) :: chunk_A, chunk_B, chunk_C

  call spamm_new_chunk(number_dimensions, N_contiguous, chunk_A)
  call spamm_new_chunk(number_dimensions, N_contiguous, chunk_B)
  call spamm_new_chunk(number_dimensions, N_contiguous, chunk_C)

  call spamm_chunk_get_N_lower(N_lower, chunk_A)
  call spamm_chunk_get_N_upper(N_upper, chunk_A)

  N_lower(1) = 128
  N_upper(1) = 256

  N_lower(2) = 512
  N_upper(2) = 768

  call spamm_chunk_get_N_lower(N_lower, chunk_B)
  call spamm_chunk_get_N_upper(N_upper, chunk_B)

  N_lower(1) = 128
  N_upper(1) = 256

  N_lower(2) = 512
  N_upper(2) = 768

  call spamm_chunk_get_matrix(A, chunk_A)
  call spamm_chunk_get_matrix(B, chunk_B)
  call spamm_chunk_get_matrix(C, chunk_C)

  do i = 1, N_contiguous
    A(i,:) = i
    B(i,:) = i
  enddo

  call spamm_chunk_print(chunk_A)

  C = matmul(A, B)

  call spamm_chunk_print(chunk_C)

end program chunk_test
