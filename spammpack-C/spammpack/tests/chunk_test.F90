program chunk_test

  use, intrinsic :: iso_c_binding

  interface

    subroutine spamm_new_chunk (ndim, N_contiguous, chunk) bind(C, name = "spamm_new_chunk_interface")
      use, intrinsic :: iso_c_binding
      integer, intent(in) :: ndim, N_contiguous
      type(c_ptr) :: chunk
    end subroutine spamm_new_chunk

    subroutine spamm_chunk_get_N_lower (cptr, chunk) bind(C, name = "spamm_chunk_get_N_lower_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: chunk, cptr
    end subroutine spamm_chunk_get_N_lower

    subroutine spamm_chunk_get_N_upper (cptr, chunk) bind(C, name = "spamm_chunk_get_N_upper_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: chunk, cptr
    end subroutine spamm_chunk_get_N_upper

    subroutine spamm_chunk_print (chunk) bind(C, name = "spamm_chunk_print_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: chunk
    end subroutine spamm_chunk_print

  end interface

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
