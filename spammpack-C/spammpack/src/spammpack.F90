module spammpack

  use spamm_derived

  IMPLICIT NONE

  interface

    subroutine spamm_new_chunk (ndim, N_contiguous, chunk) bind(C, name = "spamm_new_chunk_interface")
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: ndim, N_contiguous
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

contains

end module spammpack
