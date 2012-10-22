program chunk_test

  use :: spammpack
  use, intrinsic :: iso_c_binding

  integer :: number_dimensions = 2
  integer :: N_contiguous = 8
  integer, dimension(:), pointer :: N_lower
  integer, dimension(:), pointer :: N_upper
  real*4, dimension(:,:), pointer :: A

  type(c_ptr) :: chunk

  call spamm_new_chunk(number_dimensions, N_contiguous, chunk)

  call spamm_chunk_get_N_lower(N_lower, chunk)
  call spamm_chunk_get_N_upper(N_upper, chunk)

  N_lower(1) = 128
  N_lower(2) = 512

  N_upper(1) = 256
  N_upper(2) = 768

  call spamm_chunk_get_matrix(A, chunk)

  do i = 1, N_contiguous
    A(i,:) = i
  enddo

  call spamm_chunk_print(chunk)

end program chunk_test
