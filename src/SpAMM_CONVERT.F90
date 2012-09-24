!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 3 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SPAMMPACK)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------

!> Defines conversion operation between different data structures and SpAMM.
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
    TYPE(QuTree), POINTER                        :: qA
    INTEGER                                      :: i_lower, i_upper
    INTEGER                                      :: j_lower, j_upper

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
        !ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))

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

      ALLOCATE(qA%Quad00)
      ALLOCATE(qA%Quad01)
      ALLOCATE(qA%Quad10)
      ALLOCATE(qA%Quad11)

#ifdef OLDCONVERT
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:A_rows/2,        1:A_cols/2),        qA%Quad00)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:A_rows/2,        A_cols/2+1:A_cols), qA%Quad01)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(A_rows/2+1:A_rows, 1:A_cols/2),        qA%Quad10)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(A_rows/2+1:A_rows, A_cols/2+1:A_cols), qA%Quad11)
#else
      ! Avoid slicing here for performance.
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad00, &
        i_lower,          i_lower+A_rows/2-1, j_lower,          j_lower+A_cols/2-1)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad01, &
        i_lower,          i_lower+A_rows/2-1, j_lower+A_cols/2, j_upper)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad10, &
        i_lower+A_rows/2, i_upper,           j_lower,          j_lower+A_cols/2-1)
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A, qA%Quad11, &
        i_lower+A_rows/2, i_upper,           j_lower+A_cols/2, j_upper)
#endif

    ENDIF

  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur

END MODULE SpAMM_CONVERT
