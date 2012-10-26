program chunk_test

  use :: spammpack
  use, intrinsic :: iso_c_binding

  integer :: number_dimensions = 2
  integer :: N_contiguous = 8
  integer :: N_block = 4

  integer, dimension(:), pointer :: N
  integer, dimension(:), pointer :: N_lower
  integer, dimension(:), pointer :: N_upper

  real*4, dimension(:,:), pointer :: A, B, C
  real*4, dimension(:), pointer :: norm

  type(c_ptr) :: chunk_A, chunk_B, chunk_C

  allocate(N(number_dimensions))
  allocate(N_lower(number_dimensions))
  allocate(N_upper(number_dimensions))

  do i = 1, number_dimensions
    N(i) = N_contiguous
    N_lower(i) = 0
    N_upper(i) = N_contiguous
  enddo

  call spamm_new_chunk(number_dimensions, N_block, N, N_lower, N_upper, chunk_A)
  call spamm_new_chunk(number_dimensions, N_block, N, N_lower, N_upper, chunk_B)
  call spamm_new_chunk(number_dimensions, N_block, N, N_lower, N_upper, chunk_C)

  deallocate(N)
  deallocate(N_lower)
  deallocate(N_upper)

  call spamm_chunk_get_N_lower(N_lower, chunk_A)
  call spamm_chunk_get_N_upper(N_upper, chunk_A)

  N_lower(1) = 8
  N_upper(1) = 16

  N_lower(2) = 32
  N_upper(2) = 40

  call spamm_chunk_get_N_lower(N_lower, chunk_B)
  call spamm_chunk_get_N_upper(N_upper, chunk_B)

  N_lower(1) = 32
  N_upper(1) = 40

  N_lower(2) = 8
  N_upper(2) = 16

  call spamm_chunk_get_N_lower(N_lower, chunk_C)
  call spamm_chunk_get_N_upper(N_upper, chunk_C)

  N_lower(1) = 8
  N_upper(1) = 16

  N_lower(2) = 8
  N_upper(2) = 16

  call spamm_chunk_get_matrix(A, chunk_A)
  call spamm_chunk_get_matrix(B, chunk_B)
  call spamm_chunk_get_matrix(C, chunk_C)

  do i = 1, N_contiguous
    A(i,:) = i
    B(i,:) = i
  enddo

  C = matmul(A, B)

  call spamm_chunk_print(chunk_A)
  call spamm_chunk_print(chunk_B)
  call spamm_chunk_print(chunk_C)

end program chunk_test
