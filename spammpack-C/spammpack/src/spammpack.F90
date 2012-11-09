!> @file

#include "config_fortran.h"

module spammpack

  use, intrinsic :: iso_c_binding
  use spamm_derived

  implicit none

  interface

    !> Interface for spamm_new_chunk().
    subroutine spamm_new_chunk (ndim, N_block, N, N_lower, N_upper, chunk) &
        bind(C, name = "spamm_new_chunk_interface")
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: ndim, N_block
      integer(c_int), dimension(ndim), intent(in) :: N, N_lower, N_upper
      type(c_ptr), intent(out) :: chunk
    end subroutine spamm_new_chunk

    !> Interface for spamm_delete_chunk().
    subroutine spamm_delete_chunk (chunk) &
        bind(C, name = "spamm_delete_chunk_interface")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(inout) :: chunk
    end subroutine spamm_delete_chunk

    !> Interface for spamm_chunk_get_number_dimensions().
    subroutine spamm_chunk_get_number_dimensions (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_number_dimensions_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_number_dimensions

    !> Interface for spamm_chunk_get_N_contiguous().
    subroutine spamm_chunk_get_N_contiguous (N_contiguous, chunk) &
        bind(C, name = "spamm_chunk_get_N_contiguous_interface")
      use, intrinsic :: iso_C_binding
      integer(c_int) :: N_contiguous
      type(c_ptr) :: chunk
    end subroutine spamm_chunk_get_N_contiguous

    !> Interface for spamm_chunk_get_number_tiers().
    subroutine spamm_chunk_get_number_tiers (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_number_tiers")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_number_tiers

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

    !> Interface for spamm_chunk_get_norm().
    subroutine spamm_chunk_get_norm_interface (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_norm_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_norm_interface

    !> Interface for spamm_chunk_get_norm2().
    subroutine spamm_chunk_get_norm2_interface (cptr, chunk) &
        bind(C, name = "spamm_chunk_get_norm2_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: cptr, chunk
    end subroutine spamm_chunk_get_norm2_interface

    !> Interface for spamm_chunk_print().
    subroutine spamm_chunk_print (chunk) &
        bind(C, name = "spamm_chunk_print_interface")
      use, intrinsic :: iso_C_binding
      type(c_ptr) :: chunk
    end subroutine spamm_chunk_print

  end interface

contains

  !> Get the N_lower vector of a SpAMM chunk.
  !>
  !> @param N_lower The N_lower vector.
  !> @param chunk The chunk.
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

  !> Get the N_upper vector of a SpAMM chunk.
  !>
  !> @param N_upper The N_upper vector.
  !> @param chunk The chunk.
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

  !> Get the matrix of a SpAMM chunk.
  !>
  !> @param A The matrix.
  !> @param chunk The chunk.
  subroutine spamm_chunk_get_matrix (A, chunk)

    real*4, dimension(:,:), pointer, intent(out) :: A
    type(c_ptr), intent(in) :: chunk

    integer :: N_contiguous
    type(c_ptr) :: cptr

    call spamm_chunk_get_N_contiguous(N_contiguous, chunk)
    call spamm_chunk_get_matrix_interface(cptr, chunk)
    call c_f_pointer(cptr, A, [N_contiguous, N_contiguous])

  end subroutine spamm_chunk_get_matrix

  !> Get the norm vector of a SpAMM chunk.
  !>
  !> @param norm The norm vector.
  !> @param chunk The chunk.
  subroutine spamm_chunk_get_norm (norm, chunk)

    real*4, dimension(:), pointer, intent(out) :: norm
    type(c_ptr), intent(in) :: chunk

    integer :: N_contiguous
    type(c_ptr) :: cptr

    integer :: number_entries, N

    call spamm_chunk_get_N_contiguous(N_contiguous, chunk)
    call spamm_chunk_get_norm_interface(cptr, chunk)

    ! Figure out the number of norm entries.
    number_entries = 0
    N = N_contiguous
    do while (.true.)
      number_entries = number_entries + N_contiguous/N
      N = N/2
      if (N < SPAMM_N_BLOCK) then
        exit
      endif
    enddo

    call c_f_pointer(cptr, norm, [number_entries])

  end subroutine spamm_chunk_get_norm

  !> Get the square of the norm vector of a SpAMM chunk.
  !>
  !> @param norm2 The norm2 vector.
  !> @param chunk The chunk.
  subroutine spamm_chunk_get_norm2 (norm2, chunk)

    real*4, dimension(:), pointer, intent(out) :: norm2
    type(c_ptr), intent(in) :: chunk

    integer :: N_contiguous
    type(c_ptr) :: cptr

    integer :: number_entries, N

    call spamm_chunk_get_N_contiguous(N_contiguous, chunk)
    call spamm_chunk_get_norm2_interface(cptr, chunk)

    ! Figure out the number of norm entries.
    number_entries = 0
    N = N_contiguous
    do while (.true.)
      number_entries = number_entries + N_contiguous/N
      N = N/2
      if (N < SPAMM_N_BLOCK) then
        exit
      endif
    enddo

    call c_f_pointer(cptr, norm2, [number_entries])

  end subroutine spamm_chunk_get_norm2

end module spammpack
