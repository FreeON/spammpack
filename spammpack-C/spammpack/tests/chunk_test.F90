program chunk_test

  INTEGER :: number_dimensions = 2
  INTEGER :: N_contiguous = 128
  INTEGER*4, DIMENSION(:), POINTER :: N_lower
  INTEGER*4, DIMENSION(:), POINTER :: N_upper

  INTEGER*8 :: chunk

  call spamm_new_chunk_interface(number_dimensions, N_contiguous, chunk)

  allocate(N_lower(number_dimensions))
  allocate(N_upper(number_dimensions))

  N_lower(1) = 128
  N_lower(2) = 128

  N_upper(1) = 256
  N_upper(2) = 256

  call spamm_chunk_set_N_lower_interface(N_lower, chunk)
  call spamm_chunk_set_N_upper_interface(N_upper, chunk)

end program chunk_test
