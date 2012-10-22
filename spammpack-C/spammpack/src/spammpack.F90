!> @file

module spammpack

  use spamm_derived

  IMPLICIT NONE

  interface

    !> Interface for spamm_new_chunk().
    subroutine spamm_new_chunk (ndim, N_contiguous, chunk) bind(C, name = "spamm_new_chunk_interface")
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: ndim, N_contiguous
      type(c_ptr) :: chunk
    end subroutine spamm_new_chunk

    !> Interface for spamm_chunk_get_N_lower().
    subroutine spamm_chunk_get_N_lower (cptr, chunk) bind(C, name = "spamm_chunk_get_N_lower_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_N_lower

    !> Interface for spamm_chunk_get_N_upper().
    subroutine spamm_chunk_get_N_upper (cptr, chunk) bind(C, name = "spamm_chunk_get_N_upper_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_N_upper

    !> Interface for spamm_chunk_get_matrix().
    subroutine spamm_chunk_get_matrix (cptr, chunk) bind(C, name = "spamm_chunk_get_matrix_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_matrix

    !> Interface for spamm_chunk_get_matrix_dilated().
    subroutine spamm_chunk_get_matrix_dilated (cptr, chunk) bind(C, name = "spamm_chunk_get_matrix_dilated_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_matrix_dilated

    !> Interface for spamm_chunk_print().
    subroutine spamm_chunk_print (chunk) bind(C, name = "spamm_chunk_print_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: chunk
    end subroutine spamm_chunk_print

  end interface

contains

end module spammpack
