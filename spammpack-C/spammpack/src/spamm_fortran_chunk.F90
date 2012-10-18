MODULE SpAMM_Chunk

  IMPLICIT NONE

CONTAINS

  INTEGER*4 FUNCTION SpAMM_Chunk_get_ndim (chunk) RESULT (ndim)

    INTEGER*4, DIMENSION(1) :: chunk

    ndim = chunk(1)

  END FUNCTION SpAMM_Chunk_get_ndim

  INTEGER*4 FUNCTION SpAMM_Chunk_get_N_contiguous (chunk) RESULT (N_contiguous)

    INTEGER*4, DIMENSION(2) :: chunk

    N_contiguous = chunk(2)

  END FUNCTION SpAMM_Chunk_get_N_contiguous

  !FUNCTION SpAMM_Chunk_get_N_lower (chunk) RESULT (N_lower)

  !  INTEGER*4, DIMENSION(:) :: chunk
  !  INTEGER*4, DIMENSION(*) :: N_lower

  !END FUNCTION SpAMM_Chunk_get_N_lower

END MODULE SpAMM_Chunk
