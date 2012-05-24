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
MODULE SpAMM_MNGMENT

  USE  SpAMM_DERIVED
  USE  SpAMM_GLOBALS

  IMPLICIT NONE

  !> Interface for deep copies of SpAMM objects.
  INTERFACE Copy
    MODULE PROCEDURE SpAMM_Copy_QuTree_2_QuTree
    MODULE PROCEDURE SpAMM_Copy_QuTree_2_BiTree
    MODULE PROCEDURE SpAMM_Copy_BiTree_2_BiTree
  END INTERFACE

  !> Interface for deletion (deallocation) of SpAMM objects.
  INTERFACE Delete
    MODULE PROCEDURE SpAMM_Delete_QuTree
    MODULE PROCEDURE SpAMM_Delete_BiTree
  END INTERFACE

  !> Interface for creation (allocation) of SpAMM objects.
  INTERFACE New
    MODULE PROCEDURE SpAMM_Allocate_Full_QuTree
    MODULE PROCEDURE SpAMM_Allocate_Full_BiTree
  END INTERFACE

CONTAINS

  !=================================================================
  ! SPAMM CONTAINERS FOR MEMORY MANEGEMENT ON QUAD TREE MATRICES
  !=================================================================

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Copy QuTree into another QuTree: C <- A
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Copy_QuTree_2_QuTree(qA,qC)

    TYPE(QuTree), POINTER :: qA,qC
    INTEGER               :: Depth

    CALL NewQuNode(qC,init=.TRUE.)
    Depth=0
    !$OMP TASK UNTIED SHARED(qA,qC)
    CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA,qC,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Copy_QuTree_2_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Copy a coloumn of A into vector C: C <- Col_I(A)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Copy_QuTree_2_BiTree(qA,Col,bC)

    TYPE(QuTree), POINTER :: qA
    TYPE(BiTree), POINTER :: bC
    INTEGER               :: Col,Depth
    INTEGER,DIMENSION(2)  :: Cols

    IF(.NOT.ASSOCIATED(bC)) &
      CALL NewBiNode(bC,init=.TRUE.)
    Depth=0
    Cols=(/1,SpAMM_PADDED_MATRIX_DIMENSION/)
    !$OMP TASK UNTIED SHARED(qA,bC)
    CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA,bC,Col,Cols,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Copy_QuTree_2_BiTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Copy QuTree into another QuTree: C <- A
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Copy_BiTree_2_BiTree(bA,bC)

    TYPE(BiTree), POINTER :: bA,bC
    INTEGER               :: Depth

    IF(.NOT.ASSOCIATED(bC))&
      CALL NewBiNode(bC,init=.TRUE.)
    Depth=0
    !$OMP TASK UNTIED SHARED(bA,bC)
    CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA,bC,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Copy_BiTree_2_BiTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Delete a QuTree: A <- NULL()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Delete_QuTree(qA)

    TYPE(QuTree),POINTER :: qA
    INTEGER              :: Depth

    IF(.NOT.ASSOCIATED(qA))RETURN
    Depth=0
    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Delete_QuTree_Recur(qA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    DEALLOCATE(qA)

  END SUBROUTINE SpAMM_Delete_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Delete a BiTree: A <- NULL()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Delete_BiTree(bA)
    TYPE(BiTree),POINTER :: bA
    INTEGER              :: Depth
    IF(.NOT.ASSOCIATED(bA))RETURN
    Depth=0
    !$OMP TASK UNTIED SHARED(bA)
    CALL SpAMM_Delete_BiTree_Recur(bA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    DEALLOCATE(bA)
  END SUBROUTINE SpAMM_Delete_BiTree

  !> @brief
  !! Create a new quadtree.
  !!
  !! @details
  !! The newly created quadtree has to be deallocated by calling Delete(). If qA
  !! is already allocated then it will be free'ed by calling Delete().
  !!
  !! @param qA [inout] A pointer to a type(QuTree) object.
  SUBROUTINE SpAMM_Allocate_Full_QuTree(qA)

    TYPE(QuTree),POINTER :: qA
    INTEGER              :: Depth

    IF(ASSOCIATED(qA))CALL SpAMM_Delete_QuTree(qA)
    CALL NewQuNode(qA,init=.TRUE.)
    Depth=0
    CALL SpAMM_Allocate_Full_QuTree_Recur(qA,Depth)

  END SUBROUTINE SpAMM_Allocate_Full_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Delete a BiTree: A <- NULL()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Allocate_Full_BiTree(bA)

    TYPE(BiTree),POINTER :: bA
    INTEGER              :: Depth

    IF(ASSOCIATED(bA))CALL SpAMM_Delete_BiTree(bA)
    CALL NewBiNode(bA,init=.TRUE.)
    Depth=0
    CALL SpAMM_Allocate_Full_BiTree_Recur(bA,Depth)

  END SUBROUTINE SpAMM_Allocate_Full_BiTree

  !=================================================================
  ! SPAMM CONTAINERS FOR MEMORY MANEGEMENT
  !=================================================================
  RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur(qA,qC,Depth)

    TYPE(QuTree), POINTER  :: qA,qC
    INTEGER                :: Depth

    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(.NOT.ASSOCIATED(qC))THEN
      CALL NewQuNode(qC)
    ENDIF
    !    qC%Siz=qA%Siz
    !    qC%Lev=qA%Lev
    !    qC%Box=qA%Box
    qC%Norm=qA%Norm
    IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(qA%Blok))THEN
      !    IF(qA%Siz==SpAMM_BLOCK_SIZE.AND.ALLOCATED(qA%Blok))THEN
      IF(.NOT.ALLOCATED(qC%Blok)) &
        ALLOCATE(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
      qC%Blok=qA%Blok
    ELSE
      IF(ASSOCIATED(qA%Quad00))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad00,qC%Quad00,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad00))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad00,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad00)
      ENDIF
      IF(ASSOCIATED(qA%Quad01))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad01,qC%Quad01,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad01))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad01,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad01)
      ENDIF
      IF(ASSOCIATED(qA%Quad10))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad10,qC%Quad10,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad10))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad10,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad10)
      ENDIF
      IF(ASSOCIATED(qA%Quad11))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad11,qC%Quad11,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad11))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad11,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad11)
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Copy_BiTree_2_BiTree_Recur(bA,bC,Depth)

    TYPE(BiTree), POINTER  :: bA,bC
    INTEGER                :: Depth

    IF(.NOT.ASSOCIATED(bA))RETURN
    IF(.NOT.ASSOCIATED(bC))THEN
      CALL NewBiNode(bC)
    ENDIF
    bC%Norm=bA%Norm
    IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(bA%Vect))THEN
      IF(.NOT.ALLOCATED(bC%Vect)) &
        ALLOCATE(bC%Vect(1:SpAMM_BLOCK_SIZE))
      bC%Vect=bA%Vect
    ELSE
      IF(ASSOCIATED(bA%Sect0))THEN
        !$OMP TASK UNTIED SHARED(bA,bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA%Sect0,bC%Sect0,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(bC%Sect0))THEN
        !$OMP TASK UNTIED SHARED(bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_BiTree_Recur(bC%Sect0,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(bC%Sect0)
      ENDIF
      IF(ASSOCIATED(bA%Sect1))THEN
        !$OMP TASK UNTIED SHARED(bA,bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA%Sect1,bC%Sect1,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(bC%Sect1))THEN
        !$OMP TASK UNTIED SHARED(bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_BiTree_Recur(bC%Sect1,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(bC%Sect1)
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_BiTree_2_BiTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_BiTree_Recur(qA,bC,Col,Cols,Depth)

    TYPE(QuTree), POINTER :: qA
    TYPE(BiTree), POINTER :: bC
    INTEGER               :: I,Col,Depth,HlfSpn
    INTEGER,DIMENSION(2)  :: Cols,Col_00_10,Col_01_11

    IF(.NOT.ASSOCIATED(qA))RETURN

    IF(.NOT.ASSOCIATED(bC)) &
      CALL NewBiNode(bC)

    IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(qA%Blok))THEN
      IF(.NOT.ALLOCATED(bC%Vect)) &
        ALLOCATE(bC%Vect(1:SpAMM_BLOCK_SIZE))
      I=Col-Cols(1)+1
      bC%Vect=qA%Blok(:,I)
      !       WRITE(*,*)'I = ',I,' Col = ',Col,' Col = ',qA%Blok(:,I)
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
      bC%Vect=SpAMM_Zero
    ELSE
      HlfSpn=(Cols(2)-Cols(1))/2
      Col_00_10=(/Cols(1),Cols(1)+HlfSpn/)
      Col_01_11=(/Cols(1)+HlfSpn+1,Cols(2)/)
      IF(Col<=Col_00_10(2))THEN
        !
        IF(ASSOCIATED(qA%Quad00))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad00,bC%Sect0,Col,Col_00_10,Depth+1)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect0))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect0,Depth+1)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect0)
        ENDIF
        !
        IF(ASSOCIATED(qA%Quad10))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad10,bC%Sect1,Col,Col_00_10,Depth+1)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect1))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect1,Depth+1)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect1)
        ENDIF
        !
      ELSE
        !
        IF(ASSOCIATED(qA%Quad01))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad01,bC%Sect0,Col,Col_01_11,Depth+1)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect0))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect0,Depth+1)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect0)
        ENDIF
        !
        IF(ASSOCIATED(qA%Quad11))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad11,bC%Sect1,Col,Col_01_11,Depth+1)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect1))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect1,Depth+1)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect1)
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_QuTree_2_BiTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Delete_QuTree_Recur(qA,Depth)

    TYPE(QuTree),POINTER  :: qA
    INTEGER :: Status,Depth

    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(ALLOCATED(qA%Blok))THEN
      !$OMP CRITICAL
      DEALLOCATE(qA%Blok,STAT=Status)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(qA%Quad00))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad00,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad00)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(qA%Quad01))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad01,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad01)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(qA%Quad10))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad10,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad10)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(qA%Quad11))THEN
       !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad11,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad11)
      !$OMP END CRITICAL
    ENDIF

  END SUBROUTINE SpAMM_Delete_QuTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Delete_BiTree_Recur(bA,Depth)

    TYPE(BiTree),POINTER  :: bA
    INTEGER               :: Status,Depth

    IF(.NOT.ASSOCIATED(bA))RETURN
    IF(ALLOCATED(bA%Vect))THEN
      !$OMP CRITICAL
      DEALLOCATE(bA%Vect,STAT=Status)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(bA%Sect0))THEN
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_BiTree_Recur(bA%Sect0,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(bA%Sect0)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(bA%Sect1))THEN
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_BiTree_Recur(bA%Sect1,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(bA%Sect1)
      !$OMP END CRITICAL
    ENDIF

  END SUBROUTINE SpAMM_Delete_BiTree_Recur

  SUBROUTINE NewQuNode(qA,init)

    LOGICAL,OPTIONAL :: init
    TYPE(QuTree), POINTER :: qA

    IF(PRESENT(init))THEN
      IF(ASSOCIATED(qA))THEN
        WRITE(*,*)'LOGIC ERROR IN NewQuNode'
        CALL SpAMM_Trap()
      ENDIF
      ALLOCATE(qA)
    ELSE
      IF(.NOT.ASSOCIATED(qA))THEN
        ALLOCATE(qA)
      ELSE
        STOP ' Logic error 2 in NewQuNode '
      ENDIF
    ENDIF
    qA%Norm=SpAMM_Zero
    NULLIFY(qA%Quad00)
    NULLIFY(qA%Quad01)
    NULLIFY(qA%Quad10)
    NULLIFY(qA%Quad11)

  END SUBROUTINE NewQuNode

  SUBROUTINE NewBiNode(bA,init)

    LOGICAL,OPTIONAL :: init
    TYPE(BiTree), POINTER :: bA

    IF(PRESENT(init))THEN
      IF(ASSOCIATED(bA))STOP 'LOGIC ERROR IN NewBiNode'
      ALLOCATE(bA)
    ELSE
      IF(.NOT.ASSOCIATED(bA))THEN
        ALLOCATE(bA)
      ENDIF
    ENDIF
    bA%Norm=SpAMM_Zero
    NULLIFY(bA%Sect0)
    NULLIFY(bA%Sect1)

  END SUBROUTINE NewBiNode

  !> @private
  !!
  !! @brief
  !! Recursive allocation of a quadtree.
  !!
  !! @param qA A pointer to a type(QuTree) object.
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Allocate_Full_QuTree_Recur(qA,Depth)

    TYPE(QuTree),POINTER        :: qA
    INTEGER                     :: Depth

    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))
      qA%Blok=SpAMM_Zero
      NULLIFY(qA%Quad00)
      NULLIFY(qA%Quad01)
      NULLIFY(qA%Quad10)
      NULLIFY(qA%Quad11)
      RETURN
    ELSE
      ALLOCATE(qA%Quad00)
      ALLOCATE(qA%Quad01)
      ALLOCATE(qA%Quad10)
      ALLOCATE(qA%Quad11)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad00,Depth+1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad01,Depth+1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad10,Depth+1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad11,Depth+1)
    ENDIF

  END SUBROUTINE SpAMM_Allocate_Full_QuTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Allocate_Full_BiTree_Recur(bA,Depth)

    TYPE(BiTree),POINTER :: bA
    INTEGER              :: Depth

    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      ALLOCATE(bA%Vect(SpAMM_BLOCK_SIZE))
      bA%Vect=SpAMM_Zero
      NULLIFY(bA%Sect0)
      NULLIFY(bA%Sect1)
      RETURN
    ELSE
      ALLOCATE(bA%Sect0)
      ALLOCATE(bA%Sect1)
      CALL SpAMM_Allocate_Full_BiTree_Recur(bA%Sect0,Depth+1)
      CALL SpAMM_Allocate_Full_BiTree_Recur(bA%Sect1,Depth+1)
    ENDIF

  END SUBROUTINE SpAMM_Allocate_Full_BiTree_Recur

END MODULE SpAMM_MNGMENT
