#include "config_fortran.h"

MODULE SpAMMPACK_MANAGEMENT

  USE SpAMMPACK_TYPES

  IMPLICIT NONE

  INTERFACE SpAMM_New

    MODULE PROCEDURE SpAMM_New_SpAMM_RNK2, SpAMM_New_QuTree

  END INTERFACE SpAMM_New

  INTERFACE SpAMM_Set

    MODULE PROCEDURE SpAMM_Set_SpAMM_RNK2

  END INTERFACE SpAMM_Set

  INTERFACE SpAMM_Get

    MODULE PROCEDURE SpAMM_Get_SpAMM_RNK2

  END INTERFACE SpAMM_Get

  INTERFACE SpAMM_SetEq

    MODULE PROCEDURE SpAMM_SetEq_SpAMM_RNK2_to_Dense, &
        SpAMM_SetEq_SpAMM_C_to_Dense, &
        SpAMM_SetEq_SpAMM_C_to_SpAMM_C

  END INTERFACE SpAMM_SetEq

  INTERFACE

    !> Interface for spamm_convert_dense_to_spamm().
    SUBROUTINE spamm_convert_dense_to_spamm (ndim, N, chunk_tier, use_linear_tree, A_dense, A) &
        BIND(C, NAME = "spamm_convert_dense_to_spamm_interface")
      USE, INTRINSIC :: iso_c_binding
      INTEGER(c_int), INTENT(IN) :: ndim
      INTEGER(c_int), DIMENSION(ndim), INTENT(IN) :: N
      INTEGER(c_int), INTENT(IN) :: chunk_tier
      INTEGER(c_int), INTENT(IN) :: use_linear_tree
      REAL(KIND = c_float), DIMENSION(*) :: A_dense
      TYPE(c_ptr), INTENT(OUT) :: A
    END subroutine spamm_convert_dense_to_spamm

  END INTERFACE

CONTAINS

  SUBROUTINE SpAMM_New_QuTree (qA)

    TYPE(QuTree), INTENT(INOUT) :: qA

  END SUBROUTINE SpAMM_New_QuTree

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

  SUBROUTINE SpAMM_Set_SpAMM_RNK2 (i, j, Aij, A)

    INTEGER, INTENT(IN) :: i, j
    REAL*4, INTENT(IN) :: Aij
    TYPE(SpAMM_RNK2), INTENT(INOUT) :: A

  END SUBROUTINE SpAMM_Set_SpAMM_RNK2

  REAL*4 FUNCTION SpAMM_Get_SpAMM_RNK2 (i, j, A) RESULT(Aij)

    INTEGER, INTENT(IN) :: i, j
    TYPE(SpAMM_RNK2), INTENT(IN) :: A

    Aij = 0

  END FUNCTION SpAMM_Get_SpAMM_RNK2

  SUBROUTINE SpAMM_SetEq_SpAMM_RNK2_to_Dense (A, ADense)

    TYPE(SpAMM_RNK2), INTENT(INOUT) :: A
    REAL*4, DIMENSION(:,:), INTENT(IN) :: ADense

    INTEGER :: i, j

    IF(A%N(1) /= SIZE(ADense, 1)) THEN
      STOP
    ENDIF

    IF(A%N(2) /= SIZE(ADense, 2)) THEN
      STOP
    ENDIF

    DO i = 1, A%N(1)
      DO j = 1, A%N(2)
        CALL SpAMM_Set(i, j, ADense(i,j), A)
      ENDDO
    ENDDO

  END SUBROUTINE SpAMM_SetEq_SpAMM_RNK2_to_Dense

  SUBROUTINE SpAMM_SetEq_SpAMM_C_to_Dense (A, ADense, chunkTier, useLinearTree)

    TYPE(c_ptr), INTENT(INOUT) :: A
    REAL*8, DIMENSION(:,:), INTENT(IN) :: ADense
    INTEGER, INTENT(IN) :: chunkTier
    LOGICAL, INTENT(IN) :: useLinearTree

    REAL*4, ALLOCATABLE, DIMENSION(:,:) :: ADenseSingle
    INTEGER :: i, j
    INTEGER :: use_linear_tree

    IF(useLinearTree) THEN
      use_linear_tree = 1
    ELSE
      use_linear_tree = 0
    ENDIF

    ALLOCATE(ADenseSingle(size(ADense, 1), size(ADense, 2)))
    DO i = 1, size(ADense, 1)
      DO j = 1, size(ADense, 2)
        ADenseSingle(i,j) = ADense(i,j)
      ENDDO
    ENDDO

    CALL spamm_convert_dense_to_spamm(2, (/ size(ADense, 1), size(ADense, 2) /), &
      chunkTier, use_linear_tree, ADenseSingle, A)

    DEALLOCATE(ADenseSingle)

  END SUBROUTINE SpAMM_SetEq_SpAMM_C_to_Dense

  SUBROUTINE SpAMM_SetEq_SpAMM_C_to_SpAMM_C (A, B)

    TYPE(c_ptr), INTENT(INOUT) :: A
    TYPE(c_ptr), INTENT(IN) :: B

  END SUBROUTINE SpAMM_SetEq_SpAMM_C_to_SpAMM_C

END MODULE SpAMMPACK_MANAGEMENT
