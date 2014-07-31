!> Management functions for quadtree objects.
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
MODULE SpAMM_MNGMENT

  USE SpAMM_DERIVED
  use spamm_types

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Copy
  PUBLIC :: SpAMM_Copy_QuTree_2_QuTree_Recur
  PUBLIC :: SpAMM_Copy_BiTree_2_BiTree_Recur
  PUBLIC :: Delete
  PUBLIC :: SpAMM_Delete_QuTree_Recur
  PUBLIC :: New
  PUBLIC :: NewQuNode

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

  !> Copy QuTree into another QuTree: @f$ C \leftarrow A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qC Pointer to matrix C.
  SUBROUTINE SpAMM_Copy_QuTree_2_QuTree(qA,qC)

    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC

    INTEGER :: Depth

    CALL NewQuNode(qC)
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

  !> Delete a quadtree matrix.
  !!
  !! @param qA Pointer to quadtree matrix.
  SUBROUTINE SpAMM_Delete_QuTree(qA)

    TYPE(QuTree), POINTER :: qA

    INTEGER :: Depth

    IF(.NOT.ASSOCIATED(qA)) RETURN

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

  !> Create a new quadtree.
  !!
  !! The newly created quadtree has to be deallocated by calling Delete(). If qA
  !! is already allocated then it will be free'ed by calling Delete().
  !!
  !! @param qA A pointer to a type(QuTree) object.
  SUBROUTINE SpAMM_Allocate_Full_QuTree(qA)

    TYPE(QuTree),POINTER :: qA
    INTEGER              :: Depth

    IF(ASSOCIATED(qA)) CALL SpAMM_Delete_QuTree(qA)
    CALL NewQuNode(qA)

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

  !> Recursive copy of quadtree.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qC Pointer to matrix C.
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur (qA, qC, Depth)

    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
    INTEGER                              :: Depth

    IF(.NOT.ASSOCIATED(qA))RETURN

    IF(.NOT.ASSOCIATED(qC))THEN
      CALL NewQuNode(qC)
    ENDIF

!    qC%Norms%FrobeniusNorm=qA%Norms%FrobeniusNorm
    qC%Norm=qA%Norm
    !IF(Depth==SpAMM_TOTAL_DEPTH.AND.ASSOCIATED(qA%Blok))THEN
    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      ! IF(qA%Siz==SpAMM_BLOCK_SIZE.AND.ALLOCATED(qA%Blok))THEN
      !IF(.NOT.ASSOCIATED(qC%Blok)) THEN
      !  ALLOCATE(qC%Blok)
      !ENDIF
      qC%Blok=qA%Blok
    ELSE
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
      IF(ASSOCIATED(qA%Quad12))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad12,qC%Quad12,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad12))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad12,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad12)
      ENDIF
      IF(ASSOCIATED(qA%Quad21))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad21,qC%Quad21,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad21))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad21,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad21)
      ENDIF
      IF(ASSOCIATED(qA%Quad22))THEN
        !$OMP TASK UNTIED SHARED(qA,qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad22,qC%Quad22,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad22))THEN
        !$OMP TASK UNTIED SHARED(qC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad22,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad22)
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
    !    bC%Norms%FrobeniusNorm=bA%Norms%FrobeniusNorm
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

    !IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(qA%Blok))THEN
    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      !IF(.NOT.ALLOCATED(bC%Vect)) &
      !  ALLOCATE(bC%Vect(1:SpAMM_BLOCK_SIZE))
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
        IF(ASSOCIATED(qA%Quad11))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad11,bC%Sect0,Col,Col_00_10,Depth+1)
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
        IF(ASSOCIATED(qA%Quad21))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad21,bC%Sect1,Col,Col_00_10,Depth+1)
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
        IF(ASSOCIATED(qA%Quad12))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad12,bC%Sect0,Col,Col_01_11,Depth+1)
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
        IF(ASSOCIATED(qA%Quad22))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad22,bC%Sect1,Col,Col_01_11,Depth+1)
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

  !> Recursive part of delete quadtree.
  !!
  !! @param qA Pointer to quadtree node.
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Delete_QuTree_Recur(qA,Depth)

    TYPE(QuTree),POINTER  :: qA
    INTEGER :: Depth
    !INTEGER :: Status

    IF(.NOT.ASSOCIATED(qA))RETURN
    !$OMP CRITICAL
    !IF(ALLOCATED(qA%Blok))THEN
    !  DEALLOCATE(qA%Blok,STAT=Status)
    !ENDIF
#ifdef _OPENMP
    CALL OMP_DESTROY_LOCK(qA%lock)
#endif
    !$OMP END CRITICAL

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

    IF(ASSOCIATED(qA%Quad12))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad12,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad12)
      !$OMP END CRITICAL
    ENDIF

    IF(ASSOCIATED(qA%Quad21))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad21,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad21)
      !$OMP END CRITICAL
    ENDIF

    IF(ASSOCIATED(qA%Quad22))THEN
       !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad22,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad22)
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

  !> Allocate a new node in the quadtree.
  !!
  !! @param qA Pointer to node.
  SUBROUTINE NewQuNode(qA)

    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA

    ! Delete node if it already exists.
    IF(ASSOCIATED(qA)) THEN
      CALL Delete(qA)
    ENDIF

    ! Allocate new node.
    ALLOCATE(qA)

    ! Initialize.
    qA%number_nonzeros = 0
    qA%number_operations = 0

#ifdef _OPENMP
    ! Initialize lock.
    CALL OMP_INIT_LOCK(qA%lock)
#endif

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
    bA%Norm=SpAMM_BIG_DBL
    NULLIFY(bA%Sect0)
    NULLIFY(bA%Sect1)

  END SUBROUTINE NewBiNode

  !> @private
  !!
  !! Recursive allocation of a quadtree.
  !!
  !! @param qA A pointer to a type(QuTree) object.
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Allocate_Full_QuTree_Recur(qA, Depth)

    TYPE(QuTree), POINTER :: qA
    INTEGER               :: Depth

    ! Allocate new node.
    CALL NewQuNode(qA)

    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))
      qA%Blok=SpAMM_Zero
      RETURN
    ELSE
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad11,Depth+1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad12,Depth+1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad21,Depth+1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad22,Depth+1)
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
