MODULE SpAMMPACK_ALGEBRA

  USE SpAMMPACK_TYPES

  IMPLICIT NONE

  INTERFACE SpAMM_Multiply

    MODULE PROCEDURE SpAMM_Multiply_SpAMM_RNK2_x_SpAMM_RNK2

  END INTERFACE SpAMM_Multiply

CONTAINS

  !> Multiply a SpAMM rank 2 matrix with another one.
  SUBROUTINE SpAMM_Multiply_SpAMM_RNK2_x_SpAMM_RNK2 (tolerance, &
      alpha, A, B, beta, C)

    REAL*4, INTENT(IN) :: tolerance
    REAL*4, INTENT(IN) :: alpha
    TYPE(SpAMM_RNK2), INTENT(IN) :: A
    TYPE(SpAMM_RNK2), INTENT(IN) :: B
    TYPE(SpAMM_RNK2), INTENT(INOUT) :: C
    REAL*4, INTENT(IN) :: beta

    IF(A%chunkTier == 0) THEN
      CALL spamm_chunk_multiply_scalar(beta, C%root%chunk)
      CALL spamm_chunk_multiply(tolerance, alpha, A%root%chunk, &
        B%root%chunk, C%root%chunk)
    ELSE
      CALL SpAMM_Multiply_SpAMM_RNK2_x_Scalar_Recur(beta, C)

      !$OMP TASK UNTIED
      CALL SpAMM_Multiply_QuTree_x_QuTree_Recursive(tolerance, alpha, &
        A%root, B%root, C%root, A%N, &
        (/ 0, 0 /), (/ A%NPadded, A%NPadded /), &
        0, A%chunkTier, A%useLinearTree)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Multiply_SpAMM_RNK2_x_SpAMM_RNK2

  !> Scalar multiply: @f$ A \leftarrow a A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param a Scalar a.
  RECURSIVE SUBROUTINE SpAMM_Multiply_SpAMM_RNK2_x_Scalar (alpha, A)

    REAL*4, INTENT(IN) :: alpha
    TYPE(SpAMM_RNK2), INTENT(INOUT) :: A

    !$OMP TASK
    CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(alpha, A%root, 0, A%chunkTier, A%useLinearTree)
    !$OMP END TASK

    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Multiply_SpAMM_RNK2_x_Scalar

  !> Recursive part of scalar multiply with quadtree matrix, @f$ A \leftarrow a
  !! A @f$.
  !!
  !! @param qA Pointer to quadtree.
  !! @param a The scalar
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recursive (alpha, qA, &
      tier, chunkTier, useLinearTree)

    REAL*4, INTENT(IN) :: alpha
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    INTEGER, INTENT(IN) :: tier
    INTEGER, INTENT(IN) :: chunkTier
    LOGICAL, INTENT(IN) :: useLinearTree

    INTEGER :: i, j

    IF(.NOT.ASSOCIATED(qA)) RETURN

    IF(tier == chunkTier) THEN
      ! At the bottom, multiply the block.
      CALL spamm_chunk_multiply_scalar(alpha, qA%chunk, qA%norm2)
      qA%Norm = sqrt(qA%norm2)
    ELSE

      DO i = 1, 2
        DO j = 1, 2
          !$OMP TASK UNTIED SHARED(qA)
          CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(alpha, qA%child(i,j), tier+1, chunkTier, useLinearTree)
          !$OMP END TASK
        ENDDO
      ENDDO

      !$OMP TASKWAIT

      qA%norm2 = alpha*alpha*qA%norm2
      qA%norm = sqrt(qA%norm2)
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recursive

  !> Recursive part of the multiplication operation between a quadtree and a
  !! quadtree, @f$ C \leftarrow A \times B @f$.
  !!
  !! @param qA Pointer to quadtree A.
  !! @param qB Pointer to quadtree B.
  !! @param qC Pointer to quadtree C.
  !! @param threshold The SpAMM product tolerance.
  !! @param Depth The current Depth.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recursive (tolerance, &
      alpha, qA, qB, qC, N, NLower, NUpper, tier, chunkTier, useLinearTree)

    REAL*4, INTENT(IN) :: tolerance
    REAL*4, INTENT(IN) :: alpha
    TYPE(QuTree), POINTER, INTENT(IN) :: qA, qB
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
    INTEGER, DIMENSION(2), INTENT(IN) :: N
    INTEGER, DIMENSION(2), INTENT(IN) :: NLower
    INTEGER, DIMENSION(2), INTENT(IN) :: NUpper
    INTEGER, INTENT(IN) :: tier
    INTEGER, INTENT(IN) :: chunkTier
    LOGICAL, INTENT(IN) :: useLinearTree

    INTEGER :: i, j, k

    IF(tier == chunkTier) THEN
#ifdef _OPENMP
      CALL OMP_SET_LOCK(qC%lock)
#endif
      CALL spamm_chunk_multiply(tolerance, alpha, qA%chunk, qB%chunk, qC%chunk)
      qC%norm = sqrt(qC%norm2)
#ifdef _OPENMP
      CALL OMP_UNSET_LOCK(qC%lock)
#endif
    ELSE

      DO i = 1, 2
        DO j = 1, 2
          DO k = 1, 2

            IF(ASSOCIATED(qA%child(i,k)%node) .AND. ASSOCIATED(qB%child(k,j)%node)) THEN
              IF(qA%child(i,k)%node%norm*qB%child(k,j)%node%norm > tolerance) THEN

#ifdef _OPENMP
                CALL OMP_SET_LOCK(qC%lock)
#endif
                IF(.NOT. ASSOCIATED(qC%child(i,j)%node)) THEN
                  CALL SpAMM_New_QuTree(qC%child(i,j)%node)
                ENDIF
#ifdef _OPENMP
                CALL OMP_UNSET_LOCK(qC%lock)
#endif

                !$OMP TASK UNTIED SHARED(qA,qB,qC)
                CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(tolerance, alpha, &
                  qA%child(i,k)%node, qB%child(k,j)%node, qC%child(i,j)%node, &
                  N, NLower, NUpper, tier+1, chunkTier, useLinearTree)
                !$OMP END TASK

              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !$OMP TASKWAIT
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recursive

END MODULE SpAMMPACK_ALGEBRA
