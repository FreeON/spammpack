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

    qA=>NULL()
    CALL NewQuNode(qA, init = .TRUE.)
    CALL SpAMM_Convert_Dense_2_QuTree_Recur(A,qA)

    ! Update norms.
    qA%Norms=Norm(qA)
    qA%Norms%FrobeniusNorm = SQRT(qA%Norms%FrobeniusNorm)

  END FUNCTION SpAMM_Convert_Dense_2_QuTree

  !> Recursively convert a dense matrix to a quadtree.
  !!
  !! @param A The dense matrix.
  !! @param qA A pointer to a quadtree node.
  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur (A, qA)

    REAL(SpAMM_KIND), DIMENSION(:,:), INTENT(IN) :: A
    TYPE(QuTree), POINTER                        :: qA

    INTEGER          :: I,J
    TYPE(SpAMM_Norm) :: Norms

    I=SIZE(A,1)
    J=SIZE(A,2)
    IF(I<=SpAMM_BLOCK_SIZE.AND.J<=SpAMM_BLOCK_SIZE)THEN
      IF(I < SpAMM_BLOCK_SIZE .OR. J < SpAMM_BLOCK_SIZE) THEN
        WRITE(*, *) "LOGIC ERROR IN SpAMM: padding error"
        WRITE(*, *) "SIZE(A, 1) = ", I
        WRITE(*, *) "SIZE(A, 2) = ", J
        WRITE(*, *) "SpAMM_BLOCK_SIZE = ", SpAMM_BLOCK_SIZE
        CALL SpAMM_Exit(1)
      ELSE
        ! qA%Siz=SpAMM_BLOCK_SIZE
        ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))
        qA%Blok(1:I,1:J)=A(1:I,1:J)
        NULLIFY(qA%Quad00)
        NULLIFY(qA%Quad01)
        NULLIFY(qA%Quad10)
        NULLIFY(qA%Quad11)
      ENDIF
      RETURN
    ELSE

      ALLOCATE(qA%Quad00)
      ALLOCATE(qA%Quad01)
      ALLOCATE(qA%Quad10)
      ALLOCATE(qA%Quad11)

      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:I/2  ,1:J/2  )  , qA%Quad00 )
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:I/2  ,J/2+1:J)  , qA%Quad01 )
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(I/2+1:I,1:J/2  )  , qA%Quad10 )
      CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(I/2+1:I,J/2+1:J)  , qA%Quad11 )

    ENDIF

  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur

END MODULE SpAMM_CONVERT
