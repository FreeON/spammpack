PROGRAM create_SpAMM

  USE SpAMMPACK

  IMPLICIT NONE

  TYPE(SpAMM_RNK2) :: A
  INTEGER, DIMENSION(2) :: N = (/ 513, 513 /)
  INTEGER :: chunkTier = 2
  LOGICAL :: useLinearTree = .FALSE.

  REAL*4, ALLOCATABLE, DIMENSION(:,:) :: ADense
  INTEGER :: i, j

  CALL SpAMM_New(N, chunkTier, useLinearTree, A)

  ALLOCATE(ADense(N(1),N(2)))

  DO i = 1, N(1)
    DO j = 1, N(2)
      CALL random_number(ADense(i,j))
    ENDDO
  ENDDO

  CALL SpAMM_SetEq(A, ADense)

  DO i = 1, N(1)
    DO j = 1, N(2)
      IF(SpAMM_Get(i, j, A) /= ADense(i,j)) THEN
        WRITE(*,*) "error, found ", SpAMM_Get(i, j, A), " but should have found ", ADense(i,j)
        STOP
      ENDIF
    ENDDO
  ENDDO

END PROGRAM create_SpAMM
