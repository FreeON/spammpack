#include "config_fortran.h"

MODULE SpAMMPACK_MANAGEMENT

  USE SpAMMPACK_TYPES

  IMPLICIT NONE

  INTERFACE SpAMM_New

    MODULE PROCEDURE SpAMM_New_SpAMM_RNK2

  END INTERFACE SpAMM_New

CONTAINS

  INTEGER FUNCTION get_tree_depth (numberDimensions, N, useLinearTree) RESULT(depth)

    INTEGER, INTENT(IN) :: numberDimensions
    INTEGER, DIMENSION(numberDimensions), INTENT(IN) :: N
    LOGICAL, INTENT(IN) :: useLinearTree

    INTEGER :: NTemp, dim
    REAL*4 :: x, xN

    x = 0E0;
    DO dim = 1, numberDimensions
      IF(useLinearTree) THEN
        IF(N(dim) < SPAMM_N_KERNEL) THEN
          NTemp = SPAMM_N_KERNEL
        ELSE
          NTemp = N(dim)
        ENDIF
      ELSE
        NTemp = N(dim)
      ENDIF

      xN = log(REAL(NTemp))/log(2E0)

      IF(xN > x) THEN
        x = xN
      ENDIF
    ENDDO

    depth = ceiling(x)

    ! Double check the depth.
    IF(depth >= 1) THEN
      DO dim = 1, numberDimensions
        IF(2**(depth-1) < N(dim)) THEN
          depth = depth+1
          EXIT
        ENDIF
      ENDDO
      depth = depth-1
    ENDIF

  END FUNCTION get_tree_depth

  SUBROUTINE SpAMM_New_SpAMM_RNK2 (N, chunkTier, useLinearTree, A)

    INTEGER, DIMENSION(2), INTENT(IN) :: N
    INTEGER, INTENT(IN) :: chunkTier
    LOGICAL, INTENT(IN) :: useLinearTree
    TYPE(SpAMM_RNK2), INTENT(INOUT) :: A

    A%N = N
    A%depth = get_tree_depth(2, N, useLinearTree)
    A%NPadded = 2**A%depth
    IF(chunkTier > A%depth) THEN
      A%chunkTier = A%depth
    ELSE
      A%chunkTier = chunkTier
      A%depth = chunkTier
    ENDIF

  END SUBROUTINE SpAMM_New_SpAMM_RNK2

END MODULE SpAMMPACK_MANAGEMENT
