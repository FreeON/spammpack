!> @file

#include "config_fortran.h"

!> Chunk interface module.
MODULE SpAMMPACK_CHUNK

  USE, INTRINSIC :: iso_C_binding

  IMPLICIT NONE

  INTERFACE

    !> Interface for spamm_chunk_multiply().
    SUBROUTINE spamm_chunk_multiply (tolerance, alpha, chunkA, chunkB, chunkC, norm2) &
        BIND(C, name = "spamm_chunk_multiply_interface")
      USE, INTRINSIC :: iso_C_binding
      REAL(c_float), INTENT(IN) :: tolerance
      REAL(c_float), INTENT(IN) :: alpha
      TYPE(c_ptr), INTENT(IN) :: chunkA
      TYPE(c_ptr), INTENT(IN) :: chunkB
      TYPE(c_ptr), INTENT(INOUT) :: chunkC
      REAL(c_float), INTENT(OUT) :: norm2
    END SUBROUTINE spamm_chunk_multiply

    SUBROUTINE spamm_chunk_multiply_scalar (alpha, chunk, norm2) &
        BIND(C, name = "spamm_chunk_multiply_scalar_interface")
      USE, INTRINSIC :: iso_C_binding
      REAL(c_float), INTENT(IN) :: alpha
      TYPE(c_ptr), INTENT(INOUT) :: chunk
      REAL(c_float), INTENT(OUT) :: norm2
    END SUBROUTINE spamm_chunk_multiply_scalar

    !> Interface for spamm_new_chunk().
    SUBROUTINE spamm_new_chunk (ndim, N_block, N, N_lower, N_upper, chunk) &
        BIND(C, name = "spamm_new_chunk_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int), INTENT(IN) :: ndim, N_block
      INTEGER(c_int), DIMENSION(ndim), INTENT(IN) :: N, N_lower, N_upper
      TYPE(c_ptr), INTENT(OUT) :: chunk
    END SUBROUTINE spamm_new_chunk

    !> Interface for spamm_delete_chunk().
    SUBROUTINE spamm_delete_chunk (chunk) &
        BIND(C, name = "spamm_delete_chunk_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr), INTENT(INOUT) :: chunk
    END SUBROUTINE spamm_delete_chunk

    !> Interface for spamm_chunk_get_number_dimensions().
    SUBROUTINE spamm_chunk_get_number_dimensions (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_number_dimensions_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_number_dimensions

    !> Interface for spamm_chunk_get_N_contiguous().
    SUBROUTINE spamm_chunk_get_N_contiguous (N_contiguous, chunk) &
        BIND(C, name = "spamm_chunk_get_N_contiguous_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int) :: N_contiguous
      TYPE(c_ptr) :: chunk
    END SUBROUTINE spamm_chunk_get_N_contiguous

    !> Interface for spamm_chunk_get_number_tiers().
    SUBROUTINE spamm_chunk_get_number_tiers (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_number_tiers")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_number_tiers

    !> Interface for spamm_chunk_get_N_lower().
    SUBROUTINE spamm_chunk_get_N_lower_interface (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_N_lower_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_N_lower_interface

    !> Interface for spamm_chunk_get_N_upper().
    SUBROUTINE spamm_chunk_get_N_upper_interface (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_N_upper_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_N_upper_interface

    !> Interface for spamm_chunk_get_matrix().
    SUBROUTINE spamm_chunk_get_matrix_interface (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_matrix_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_matrix_interface

    !> Interface for spamm_chunk_get_matrix_dilated().
    SUBROUTINE spamm_chunk_get_matrix_dilated (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_matrix_dilated_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_matrix_dilated

    !> Interface for spamm_chunk_get_norm().
    SUBROUTINE spamm_chunk_get_norm_interface (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_norm_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_norm_interface

    !> Interface for spamm_chunk_get_norm2().
    SUBROUTINE spamm_chunk_get_norm2_interface (cptr, chunk) &
        BIND(C, name = "spamm_chunk_get_norm2_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: cptr, chunk
    END SUBROUTINE spamm_chunk_get_norm2_interface

    !> Interface for spamm_print_chunk().
    SUBROUTINE spamm_print_chunk (chunk) &
        BIND(C, name = "spamm_print_chunk_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr) :: chunk
    END SUBROUTINE spamm_print_chunk

  END INTERFACE

CONTAINS

  !> Get the N_lower vector of a SpAMM chunk.
  !>
  !> @param N_lower The N_lower vector.
  !> @param chunk The chunk.
  SUBROUTINE spamm_chunk_get_N_lower (N_lower, chunk)

    INTEGER(c_int), DIMENSION(:), POINTER, INTENT(OUT) :: N_lower
    TYPE(c_ptr), INTENT(IN) :: chunk

    INTEGER(c_int), POINTER :: number_dimensions
    TYPE(c_ptr) :: cptr

    CALL spamm_chunk_get_number_dimensions(cptr, chunk)
    CALL c_f_pointer(cptr, number_dimensions)

    CALL spamm_chunk_get_N_lower_interface(cptr, chunk)
    CALL c_f_pointer(cptr, N_lower, [number_dimensions])

  END SUBROUTINE spamm_chunk_get_N_lower

  !> Get the N_upper vector of a SpAMM chunk.
  !>
  !> @param N_upper The N_upper vector.
  !> @param chunk The chunk.
  SUBROUTINE spamm_chunk_get_N_upper (N_upper, chunk)

    INTEGER(c_int), DIMENSION(:), POINTER, INTENT(OUT) :: N_upper
    TYPE(c_ptr), INTENT(IN) :: chunk

    INTEGER(c_int), POINTER :: number_dimensions
    TYPE(c_ptr) :: cptr

    CALL spamm_chunk_get_number_dimensions(cptr, chunk)
    CALL c_f_pointer(cptr, number_dimensions)

    CALL spamm_chunk_get_N_upper_interface(cptr, chunk)
    CALL c_f_pointer(cptr, N_upper, [number_dimensions])

  END SUBROUTINE spamm_chunk_get_N_upper

  !> Get the matrix of a SpAMM chunk.
  !>
  !> @param A The matrix.
  !> @param chunk The chunk.
  SUBROUTINE spamm_chunk_get_matrix (A, chunk)

    REAL*4, DIMENSION(:,:), POINTER, INTENT(OUT) :: A
    TYPE(c_ptr), INTENT(IN) :: chunk

    INTEGER(c_int) :: N_contiguous
    TYPE(c_ptr) :: cptr

    CALL spamm_chunk_get_N_contiguous(N_contiguous, chunk)
    CALL spamm_chunk_get_matrix_interface(cptr, chunk)
    CALL c_f_pointer(cptr, A, [N_contiguous, N_contiguous])

  END SUBROUTINE spamm_chunk_get_matrix

  !> Get the norm vector of a SpAMM chunk.
  !>
  !> @param norm The norm vector.
  !> @param chunk The chunk.
  SUBROUTINE spamm_chunk_get_norm (norm, chunk)

    REAL*4, DIMENSION(:), POINTER, INTENT(OUT) :: norm
    TYPE(c_ptr), INTENT(IN) :: chunk

    INTEGER(c_int) :: N_contiguous
    TYPE(c_ptr) :: cptr

    INTEGER(c_int) :: number_entries, N

    CALL spamm_chunk_get_N_contiguous(N_contiguous, chunk)
    CALL spamm_chunk_get_norm_interface(cptr, chunk)

    ! Figure out the number of norm entries.
    number_entries = 0
    N = N_contiguous
    DO WHILE (.true.)
      number_entries = number_entries + N_contiguous/N
      N = N/2
      IF (N < SPAMM_N_BLOCK) THEN
        exit
      ENDIF
    ENDDO

    CALL c_f_pointer(cptr, norm, [number_entries])

  END SUBROUTINE spamm_chunk_get_norm

  !> Get the square of the norm vector of a SpAMM chunk.
  !>
  !> @param norm2 The norm2 vector.
  !> @param chunk The chunk.
  SUBROUTINE spamm_chunk_get_norm2 (norm2, chunk)

    REAL(c_float), DIMENSION(:), POINTER, INTENT(OUT) :: norm2
    TYPE(c_ptr), INTENT(IN) :: chunk

    INTEGER(c_int) :: N_contiguous
    TYPE(c_ptr) :: cptr

    INTEGER(c_int) :: number_entries, N

    CALL spamm_chunk_get_N_contiguous(N_contiguous, chunk)
    CALL spamm_chunk_get_norm2_interface(cptr, chunk)

    ! Figure out the number of norm entries.
    number_entries = 0
    N = N_contiguous
    DO WHILE (.true.)
      number_entries = number_entries + N_contiguous/N
      N = N/2
      IF (N < SPAMM_N_BLOCK) THEN
        exit
      ENDIF
    ENDDO

    CALL c_f_pointer(cptr, norm2, [number_entries])

  END SUBROUTINE spamm_chunk_get_norm2

END MODULE SpAMMPACK_CHUNK
