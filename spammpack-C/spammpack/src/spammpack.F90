!> @file

module spammpack

  use, intrinsic :: iso_c_binding
  use spamm_derived

  implicit none

  interface

    !> Interface for spamm_new_chunk().
    subroutine spamm_new_chunk (ndim, N_contiguous, chunk) &
        bind(C, name = "spamm_new_chunk_interface")
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: ndim, N_contiguous
      type(c_ptr), intent(out) :: chunk
    end subroutine spamm_new_chunk

    !> Interface for spamm_chunk_get_number_dimensions().
    subroutine spamm_chunk_get_number_dimensions (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_number_dimensions_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_number_dimensions

    !> Interface for spamm_chunk_get_N_contiguous().
    subroutine spamm_chunk_get_N_contiguous (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_N_contiguous_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_N_contiguous

    !> Interface for spamm_chunk_get_N_lower().
    subroutine spamm_chunk_get_N_lower_interface (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_N_lower_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_N_lower_interface

    !> Interface for spamm_chunk_get_N_upper().
    subroutine spamm_chunk_get_N_upper_interface (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_N_upper_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_N_upper_interface

    !> Interface for spamm_chunk_get_matrix().
    subroutine spamm_chunk_get_matrix_interface (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_matrix_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_matrix_interface

    !> Interface for spamm_chunk_get_matrix_dilated().
    subroutine spamm_chunk_get_matrix_dilated (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_matrix_dilated_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_matrix_dilated

    !> Interface for spamm_chunk_print().
    subroutine spamm_chunk_print (chunk) &
        bind(C, name = "spamm_chunk_print_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: chunk
    end subroutine spamm_chunk_print

  end interface

contains

  subroutine spamm_chunk_get_N_lower (N_lower, chunk)

    integer, dimension(:), pointer, intent(out) :: N_lower
    type(c_ptr), intent(in) :: chunk

    integer, pointer :: number_dimensions
    type(c_ptr) :: cptr

    call spamm_chunk_get_number_dimensions(cptr, chunk)
    call c_f_pointer(cptr, number_dimensions)

    call spamm_chunk_get_N_lower_interface(cptr, chunk)
    call c_f_pointer(cptr, N_lower, [number_dimensions])

  end subroutine spamm_chunk_get_N_lower

  subroutine spamm_chunk_get_N_upper (N_upper, chunk)

    integer, dimension(:), pointer, intent(out) :: N_upper
    type(c_ptr), intent(in) :: chunk

    integer, pointer :: number_dimensions
    type(c_ptr) :: cptr

    call spamm_chunk_get_number_dimensions(cptr, chunk)
    call c_f_pointer(cptr, number_dimensions)

    call spamm_chunk_get_N_upper_interface(cptr, chunk)
    call c_f_pointer(cptr, N_upper, [number_dimensions])

  end subroutine spamm_chunk_get_N_upper

  subroutine spamm_chunk_get_matrix (A, chunk)

    real*4, dimension(:,:), pointer, intent(out) :: A
    type(c_ptr), intent(in) :: chunk

    integer, pointer :: N_contiguous
    type(c_ptr) :: cptr

    call spamm_chunk_get_N_contiguous(cptr, chunk)
    call c_f_pointer(cptr, N_contiguous)

    call spamm_chunk_get_matrix_interface(cptr, chunk)
    call c_f_pointer(cptr, A, [N_contiguous, N_contiguous])

  end subroutine spamm_chunk_get_matrix

end module spammpack
