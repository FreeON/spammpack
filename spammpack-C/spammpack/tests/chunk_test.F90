program chunk_test

  INTEGER :: number_dimensions = 2
  INTEGER :: N_contiguous = 128
  INTEGER*4, DIMENSION(:), POINTER :: N_lower

  INTEGER*8 :: chunk

  chunk = spamm_new_chunk(number_dimensions, N_contiguous)
  N_lower = spamm_chunk_get_N_lower(chunk)

end program chunk_test
