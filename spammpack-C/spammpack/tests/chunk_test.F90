program chunk_test

  use :: spammpack
  use, intrinsic :: iso_c_binding

  INTEGER :: number_dimensions = 2
  INTEGER :: N_contiguous = 128
  INTEGER*4, DIMENSION(:), POINTER :: N_lower
  INTEGER*4, DIMENSION(:), POINTER :: N_upper

  type(c_ptr) :: chunk
  type(c_ptr) :: cptr

  call spamm_new_chunk(number_dimensions, N_contiguous, chunk)

  call spamm_chunk_get_N_lower(cptr, chunk)
  call c_f_pointer(cptr, N_lower, [number_dimensions])

  call spamm_chunk_get_N_upper(cptr, chunk)
  call c_f_pointer(cptr, N_upper, [number_dimensions])

  N_lower(1) = 128
  N_lower(2) = 512

  N_upper(1) = 256
  N_upper(2) = 768

  call spamm_chunk_print(chunk)

end program chunk_test
