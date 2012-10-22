program chunk_test

  use iso_c_binding

  interface

    function spamm_new_chunk (ndim, N_contiguous) result (chunk) bind(C, name = "spamm_new_chunk")
      use iso_c_binding
      integer, intent(in) :: ndim, N_contiguous
      type(c_ptr) :: chunk
    end function spamm_new_chunk

    function spamm_chunk_get_N_lower (chunk) result (cptr) bind(C, name = "spamm_chunk_get_N_lower")
      use iso_C_binding
      type(c_ptr) :: chunk, cptr
    end function spamm_chunk_get_N_lower

    function spamm_chunk_get_N_upper (chunk) result (cptr) bind(C, name = "spamm_chunk_get_N_upper")
      use iso_C_binding
      type(c_ptr) :: chunk, cptr
    end function spamm_chunk_get_N_upper

    subroutine spamm_chunk_print (chunk) bind(C, name = "spamm_chunk_print")
      use iso_C_binding
      type(c_ptr) :: chunk
    end subroutine spamm_chunk_print

  end interface

  INTEGER :: number_dimensions = 2
  INTEGER :: N_contiguous = 128
  INTEGER*4, DIMENSION(:), POINTER :: N_lower
  INTEGER*4, DIMENSION(:), POINTER :: N_upper

  type(c_ptr) :: chunk
  type(c_ptr) :: cptr

  chunk = spamm_new_chunk(number_dimensions, N_contiguous)

  cptr = spamm_chunk_get_N_lower(chunk)
  call c_f_pointer(cptr, N_lower, [number_dimensions])

  cptr = spamm_chunk_get_N_upper(chunk)
  call c_f_pointer(cptr, N_upper, [number_dimensions])

  N_lower(1) = 128
  N_lower(2) = 512

  N_upper(1) = 256
  N_upper(2) = 768

  call spamm_chunk_print(chunk)

end program chunk_test
