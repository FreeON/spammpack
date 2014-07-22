!> Defines conversion operation between different data structures and SpAMM.
!!
!! @copyright
!!
!! Copyright (c) 2014, Los Alamos National Laboratory
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice, this
!! list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!! this list of conditions and the following disclaimer in the documentation
!! and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its contributors
!! may be used to endorse or promote products derived from this software without
!! specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! @author Matt Challacombe matt.challacombe@freeon.org
!! @author Nicolas Bock nicolas.bock@freeon.org
MODULE SpAMM_CONVERT

  USE SpAMM_DERIVED
  USE SpAMM_GLOBALS
  USE SpAMM_ALGEBRA
  USE SpAMM_MNGMENT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SpAMM_Convert_Dense_2_QuTree

  CONTAINS

  !> Convert a dense matrix into a quadtree.
  !!
  !! @param A The dense matrix.
  !!
  !! @return Pointer to a new quadtree.
  FUNCTION SpAMM_Convert_Dense_2_QuTree (A) RESULT(qA)
    REAL(SpAMM_KIND), DIMENSION(:,:), INTENT(IN) :: A
    TYPE(QuTree), POINTER                        :: qA
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    ! Sanity check.
    IF(SIZE(A, 1) /= SpAMM_MATRIX_DIMENSION) THEN
      WRITE(*, *) "size(A, 1) ", SIZE(A, 1), " is different than size used for init ", SpAMM_MATRIX_DIMENSION
      STOP
    ENDIF

    IF(SIZE(A, 2) /= SpAMM_MATRIX_DIMENSION) THEN
      WRITE(*, *) "size(A, 2) ", SIZE(A, 2), " is different than size used for init ", SpAMM_MATRIX_DIMENSION
      STOP
    ENDIF

    TInitial = SpAMM_Get_Time()
    qA=>NULL()
    CALL NewQuNode(qA)
#ifdef OLDCONVERT
    CALL SpAMM_Convert_Dense_2_QuTree_Recur(A,qA)
#else
    CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA, &
      1, SpAMM_PADDED_MATRIX_DIMENSION, 1, SpAMM_PADDED_MATRIX_DIMENSION)
#endif
    ! Update norms.
    qA%Norm=Norm(qA)
    qA%Norm=SQRT(qA%Norm)
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal, "SpAMM_Convert_Dense_2_QuTree", 20)
  END FUNCTION SpAMM_Convert_Dense_2_QuTree

  !> Recursively convert a dense matrix to a quadtree.
  !!
  !! @param A The dense matrix.
  !! @param qA A pointer to a quadtree node.
  !! @param i_lower The lower value of the row index.
  !! @param i_upper The upper value of the row index.
  !! @param j_lower The lower value of the column index.
  !! @param j_upper The upper value of the column index.
#ifdef OLDCONVERT
  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur (A, qA)
#else
  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur (A, qA, &
      i_lower, i_upper, j_lower, j_upper)
#endif

    REAL(SpAMM_KIND), DIMENSION(:,:), INTENT(IN) :: A
    TYPE(QuTree), POINTER, INTENT(INOUT)         :: qA
    INTEGER, INTENT(IN)                          :: i_lower, i_upper
    INTEGER, INTENT(IN)                          :: j_lower, j_upper

    INTEGER :: i, j
    INTEGER :: A_rows, A_cols

#ifdef OLDCONVERT
    A_rows = SIZE(A,1)
    A_cols = SIZE(A,2)
#else
    A_rows = i_upper-i_lower+1
    A_cols = j_upper-j_lower+1
#endif

    IF(A_rows<=SpAMM_BLOCK_SIZE.AND.A_cols<=SpAMM_BLOCK_SIZE)THEN
      IF(A_rows < SpAMM_BLOCK_SIZE .OR. A_cols < SpAMM_BLOCK_SIZE) THEN
        WRITE(*, *) "LOGIC ERROR IN SpAMM: padding error"
        WRITE(*, *) "A_rows = ", A_rows
        WRITE(*, *) "A_cols = ", A_cols
        WRITE(*, *) "SpAMM_BLOCK_SIZE = ", SpAMM_BLOCK_SIZE
        CALL SpAMM_Exit(1)
      ELSE
        ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))

        ! Set new block to zero.
        qA%Blok = SpAMM_Zero

#ifdef OLDCONVERT
        ! This is wrong if the padded dimensions do not line up with the
        ! original matrix dimension.
        qA%Blok(1:A_rows, 1:A_cols) = A(1:A_rows, 1:A_cols)
#else
        DO i = i_lower, i_upper
          DO j = j_lower, j_upper
            ! We have to  be careful not to copy too much of A.
            IF(i <= SpAMM_MATRIX_DIMENSION .AND. j <= SpAMM_MATRIX_DIMENSION) THEN
              qA%Blok(i-i_lower+1, j-j_lower+1) = A(i, j)
            ENDIF
          ENDDO
        ENDDO
#endif
      ENDIF
      RETURN
    ELSE

      ALLOCATE(qA%Quad11)
      ALLOCATE(qA%Quad12)
      ALLOCATE(qA%Quad21)
      ALLOCATE(qA%Quad22)

#ifdef OLDCONVERT
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:A_rows/2,        1:A_cols/2),        qA%Quad11)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:A_rows/2,        A_cols/2+1:A_cols), qA%Quad12)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(A_rows/2+1:A_rows, 1:A_cols/2),        qA%Quad21)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(A_rows/2+1:A_rows, A_cols/2+1:A_cols), qA%Quad22)
#else
      ! Avoid slicing here for performance.
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad11, &
        i_lower,          i_lower+A_rows/2-1, j_lower,          j_lower+A_cols/2-1)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad12, &
        i_lower,          i_lower+A_rows/2-1, j_lower+A_cols/2, j_upper)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad21, &
        i_lower+A_rows/2, i_upper,            j_lower,          j_lower+A_cols/2-1)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad22, &
        i_lower+A_rows/2, i_upper,            j_lower+A_cols/2, j_upper)
#endif

    ENDIF

  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur

END MODULE SpAMM_CONVERT
