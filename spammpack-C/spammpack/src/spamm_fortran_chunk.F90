MODULE SpAMM_chunk

  IMPLICIT NONE

CONTAINS

  FUNCTION SpAMM_chunk_get_ndim (chunk) RESULT (ndim)

    INTEGER*8 :: chunk
    INTEGER*4, POINTER :: ndim

    ndim = chunk(1)

  END FUNCTION SpAMM_chunk_get_ndim

  FUNCTION SpAMM_chunk_get_N_contiguous (chunk) RESULT (N_contiguous)

    INTEGER*4, DIMENSION(2) :: chunk
    INTEGER*4, POINTER :: N_contiguous

    N_contiguous = chunk(2)

  END FUNCTION SpAMM_chunk_get_N_contiguous

  FUNCTION SpAMM_chunk_get_N_lower (chunk) RESULT (N_lower)

    INTEGER*4, DIMENSION(:) :: chunk
    INTEGER*4, POINTER, DIMENSION(:,:) :: N_lower

    N_lower = chunk(3)

  END FUNCTION SpAMM_chunk_get_N_lower

END MODULE SpAMM_chunk
