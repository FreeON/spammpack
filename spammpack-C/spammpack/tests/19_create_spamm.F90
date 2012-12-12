PROGRAM create_SpAMM

  USE SpAMMPACK

  IMPLICIT NONE

  TYPE(SpAMM_RNK2) :: A
  INTEGER, DIMENSION(2) :: N = (/ 513, 513 /)
  INTEGER :: chunkTier = 2
  LOGICAL :: useLinearTree = .FALSE.

  CALL SpAMM_New(N, chunkTier, useLinearTree, A)

END PROGRAM create_SpAMM
