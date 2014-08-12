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
module spamm_management

  use spamm_types
  use spamm_globals

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Copy
  PUBLIC :: SpAMM_Copy_QuTree_2_QuTree_Recur
  PUBLIC :: SpAMM_Copy_BiTree_2_BiTree_Recur
  PUBLIC :: Delete
  PUBLIC :: SpAMM_Delete_QuTree_Recur
  public :: spamm_allocate_matrix_2nd_order
  PUBLIC :: New
  PUBLIC :: NewQuNode
  public :: spamm_zero_matrix
  public :: get

  !> Interface for deep copies of SpAMM objects.
  INTERFACE Copy
    MODULE PROCEDURE SpAMM_Copy_QuTree_2_QuTree
    MODULE PROCEDURE SpAMM_Copy_QuTree_2_BiTree
    MODULE PROCEDURE SpAMM_Copy_BiTree_2_BiTree
    module procedure spamm_copy_2nd_order_to_2nd_order
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

  !> Interface for getting single matrix elements of SpAMM objects.
  interface get
    module procedure spamm_get_matrix_2nd_order
  end interface get

CONTAINS

  !> Copy a 2nd order matrix.
  !!
  !! @param A The matrix A.
  !! @param C The matrix C.
  subroutine spamm_copy_2nd_order_to_2nd_order (A, B)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    type(spamm_matrix_2nd_order), pointer, intent(inout) :: B

    if(.not. associated(A)) then
      return
    endif

    if(.not. associated(B)) then
      B => spamm_allocate_matrix_2nd_order(A%M, A%N)
    endif

    call spamm_copy_qutree_2_qutree(A%root, B%root)

  end subroutine spamm_copy_2nd_order_to_2nd_order

  !> Copy QuTree into another QuTree: @f$ C \leftarrow A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qC Pointer to matrix C.
  SUBROUTINE SpAMM_Copy_QuTree_2_QuTree(qA,qC)

    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC

    INTEGER :: Depth

    CALL NewQuNode(qC, qA%i_lower, qA%j_lower, qA%i_upper, qA%j_upper)

    !$OMP TASK UNTIED SHARED(qA,qC)
    CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA, qC)
    !$OMP END TASK
    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Copy_QuTree_2_QuTree

  !> Copy a coloumn of A into vector C: \f$ C \leftarrow A(:, j) \f$.
  !!
  !! @param qA Pointer to quadtree.
  !! @param j The column index.
  !! @param bC Pointer to bitree.
  SUBROUTINE SpAMM_Copy_QuTree_2_BiTree (qA, j, bC)

    TYPE(QuTree), POINTER, intent(in) :: qA
    TYPE(BiTree), POINTER, intent(inout) :: bC
    INTEGER, intent(in) :: j
    INTEGER,DIMENSION(2) :: Cols

    IF(.NOT.ASSOCIATED(bC)) then
      CALL NewBiNode(bC, init = .TRUE.)
    endif

    !$OMP TASK UNTIED SHARED(qA,bC)
    CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA, bC, j)
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

    TYPE(QuTree), POINTER, intent(inout) :: qA

    IF(.NOT.ASSOCIATED(qA)) RETURN

    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Delete_QuTree_Recur(qA)
    !$OMP END TASK

    !$OMP TASKWAIT
    DEALLOCATE(qA)

  END SUBROUTINE SpAMM_Delete_QuTree

  !> Delete a BiTree: \f$ A \leftarrow NULL() \f$.
  !!
  !! @param bA A pointer to a bitree.
  SUBROUTINE SpAMM_Delete_BiTree(bA)
    TYPE(BiTree),POINTER :: bA
    INTEGER              :: Depth
    IF(.NOT.ASSOCIATED(bA))RETURN
    Depth=0
    !$OMP TASK UNTIED SHARED(bA)
    CALL SpAMM_Delete_BiTree_Recur(bA)
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
  !! @param i_lower The lower row index.
  !! @param j_lower The lower column index.
  !! @param i_upper The upper row index.
  !! @param j_upper The upper column index.
  SUBROUTINE SpAMM_Allocate_Full_QuTree(qA, i_lower, j_lower, i_upper, j_upper)

    TYPE(QuTree), POINTER :: qA
    integer, intent(in) :: i_lower, j_lower, i_upper, j_upper

    IF(ASSOCIATED(qA)) CALL SpAMM_Delete_QuTree(qA)
    CALL SpAMM_Allocate_Full_QuTree_Recur(qA, i_lower, j_lower, i_upper, j_upper)

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
  RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur (qA, qC)

    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC

    IF(.NOT.ASSOCIATED(qA)) RETURN

    IF(.NOT.ASSOCIATED(qC))THEN
      CALL NewQuNode(qC, qA%i_lower, qA%j_lower, qA%i_upper, qA%j_upper)
    ENDIF

    ! qC%Norms%FrobeniusNorm=qA%Norms%FrobeniusNorm
    qC%Norm=qA%Norm
    !IF(Depth==SpAMM_TOTAL_DEPTH.AND.ASSOCIATED(qA%Blok))THEN
    IF(qA%i_upper-qA%i_lower == SPAMM_BLOCK_SIZE) then
      IF(.NOT. allocated(qC%Blok)) THEN
        ALLOCATE(qC%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
      ENDIF
      qC%Blok = qA%Blok
    ELSE
      IF(ASSOCIATED(qA%Quad11))THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad11, qC%Quad11)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad11))THEN
        !$OMP TASK UNTIED SHARED(qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad11)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad11)
      ENDIF
      IF(ASSOCIATED(qA%Quad12))THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad12, qC%Quad12)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad12))THEN
        !$OMP TASK UNTIED SHARED(qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad12)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad12)
      ENDIF
      IF(ASSOCIATED(qA%Quad21))THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad21, qC%Quad21)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad21))THEN
        !$OMP TASK UNTIED SHARED(qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad21)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(qC%Quad21)
      ENDIF
      IF(ASSOCIATED(qA%Quad22))THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad22, qC%Quad22)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(qC%Quad22))THEN
        !$OMP TASK UNTIED SHARED(qC)
        !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_QuTree_Recur(qC%Quad22)
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
        ALLOCATE(bC%Vect(1:SPAMM_BLOCK_SIZE))
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
        CALL SpAMM_Delete_BiTree_Recur(bC%Sect0)
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
        CALL SpAMM_Delete_BiTree_Recur(bC%Sect1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(bC%Sect1)
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_BiTree_2_BiTree_Recur

  !> Copy quadtree column to bitree.
  !!
  !! @param qA Pointer to quadtree.
  !! @param bC Pointer to bitree.
  !! @param j The column index.
  RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_BiTree_Recur(qA, bC, j)

    TYPE(QuTree), POINTER, intent(in) :: qA
    TYPE(BiTree), POINTER, intent(inout) :: bC
    integer, intent(in) :: j

    INTEGER               :: I, Depth, half
    INTEGER,DIMENSION(2)  :: Cols,Col_00_10,Col_01_11

    IF(.NOT.ASSOCIATED(qA)) RETURN

    IF(.NOT.ASSOCIATED(bC)) then
      CALL NewBiNode(bC)
    endif

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      if(allocated(qA%blok)) then
        IF(.NOT. ALLOCATED(bC%Vect)) then
          ALLOCATE(bC%Vect(SPAMM_BLOCK_SIZE))
        endif
        bC%Vect = qA%Blok(:, j-qA%i_lower+1)
      else
        bC%Vect = SpAMM_Zero
      endif
    ELSE
      half = (Cols(2)-Cols(1))/2
      Col_00_10 = (/ Cols(1), Cols(1)+half /)
      Col_01_11 = (/ Cols(1)+half+1, Cols(2) /)

      IF(j <= qA%i_lower+half) THEN
        IF(ASSOCIATED(qA%Quad11))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad11, bC%Sect0, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect0))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect0)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect0)
        ENDIF

        IF(ASSOCIATED(qA%Quad21))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad21, bC%Sect1, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect1))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect1)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect1)
        ENDIF
      ELSE
        IF(ASSOCIATED(qA%Quad12))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad12, bC%Sect0, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect0))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect0)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect0)
        ENDIF

        IF(ASSOCIATED(qA%Quad22))THEN
          !$OMP TASK UNTIED SHARED(qA,bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad22, bC%Sect1, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%Sect1))THEN
          !$OMP TASK UNTIED SHARED(bC) &
          !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
          CALL SpAMM_Delete_BiTree_Recur(bC%Sect1)
          !$OMP END TASK
          !$OMP TASKWAIT
          DEALLOCATE(bC%Sect1)
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_QuTree_2_BiTree_Recur

  !> Recursive part of delete quadtree.
  !!
  !! @bug Instaed of a critical section, use locks.
  !!
  !! @param qA Pointer to quadtree node.
  RECURSIVE SUBROUTINE SpAMM_Delete_QuTree_Recur(qA)

    TYPE(QuTree),POINTER  :: qA

    IF(.NOT.ASSOCIATED(qA))RETURN
    !$OMP CRITICAL
    IF(ALLOCATED(qA%Blok))THEN
      DEALLOCATE(qA%Blok)
    ENDIF
#ifdef _OPENMP
    CALL OMP_DESTROY_LOCK(qA%lock)
#endif
    !$OMP END CRITICAL

    IF(ASSOCIATED(qA%Quad11))THEN
      !$OMP TASK UNTIED SHARED(qA)
      !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad11)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad11)
      !$OMP END CRITICAL
    ENDIF

    IF(ASSOCIATED(qA%Quad12))THEN
      !$OMP TASK UNTIED SHARED(qA)
      !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad12)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad12)
      !$OMP END CRITICAL
    ENDIF

    IF(ASSOCIATED(qA%Quad21))THEN
      !$OMP TASK UNTIED SHARED(qA)
      !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad21)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad21)
      !$OMP END CRITICAL
    ENDIF

    IF(ASSOCIATED(qA%Quad22))THEN
      !$OMP TASK UNTIED SHARED(qA)
      !!$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA%Quad22)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(qA%Quad22)
      !$OMP END CRITICAL
    ENDIF

  END SUBROUTINE SpAMM_Delete_QuTree_Recur

  !> Delete a bitree.
  !!
  !! @bug Instead of a critical section, use locks.
  !!
  !! @param bA A pointer to a bitree.
  RECURSIVE SUBROUTINE SpAMM_Delete_BiTree_Recur(bA)

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
      CALL SpAMM_Delete_BiTree_Recur(bA%Sect0)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(bA%Sect0)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(bA%Sect1))THEN
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_BiTree_Recur(bA%Sect1)
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
  !! @param i_lower The lower row index.
  !! @param j_lower The lower column index.
  !! @param i_upper The upper row index.
  !! @param j_upper The upper column index.
  SUBROUTINE NewQuNode (qA, i_lower, j_lower, i_upper, j_upper)

    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    integer, intent(in) :: i_lower, j_lower, i_upper, j_upper

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
  !! @param i_lower The lower row index.
  !! @param j_lower The lower column index.
  !! @param i_upper The upper row index.
  !! @param j_upper The upper column index.
  RECURSIVE SUBROUTINE SpAMM_Allocate_Full_QuTree_Recur(qA, i_lower, j_lower, i_upper, j_upper)

    TYPE(QuTree), POINTER :: qA
    integer, intent(in) :: i_lower, j_lower, i_upper, j_upper
    integer :: rows, columns

    rows = i_upper-i_lower+1
    columns = j_upper-j_lower+1

    if(rows /= columns) then
      write(*, *) "non-square submatrix"
      error stop
    endif

    if(rows < SPAMM_BLOCK_SIZE .or. columns < SPAMM_BLOCK_SIZE) then
      write(*, *) "logic error"
      error stop
    endif

    ! Allocate new node.
    CALL NewQuNode(qA, i_lower, j_lower, i_upper, j_upper)

    IF(rows == SPAMM_BLOCK_SIZE) then
      ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
      qA%Blok = SpAMM_Zero
      RETURN
    ELSE
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad11, i_lower, j_lower, i_lower+rows/2, j_lower+columns/2)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad12, i_lower, j_lower+columns/2, i_lower+rows/2, j_upper)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad21, i_lower+rows/2, j_lower, i_upper, j_lower+columns/2)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad22, i_lower+rows/2, j_lower+columns/2, i_upper, j_upper)
    ENDIF

  END SUBROUTINE SpAMM_Allocate_Full_QuTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Allocate_Full_BiTree_Recur(bA,Depth)

    TYPE(BiTree),POINTER :: bA
    INTEGER              :: Depth

    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      ALLOCATE(bA%Vect(SPAMM_BLOCK_SIZE))
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

  !> Construct a zero matrix of size M x N.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !!
  !! @return The matrix.
  function spamm_zero_matrix (M, N) result (A)

    type(spamm_matrix_2nd_order), pointer :: A
    integer, intent(in) :: M, N

    A => spamm_allocate_matrix_2nd_order(M, N)
    call spamm_allocate_full_qutree(A%root, 1, 1, A%N_padded, A%N_padded)

  end function spamm_zero_matrix

  !> Get a matrix element from a 2nd order SpAMM matrix.
  !!
  !! @bug This function is not implemented yet.
  !!
  !! @param A The matrix.
  !! @param i The row index.
  !! @param j The column index.
  !!
  !! @return The matrix element.
  function spamm_get_matrix_2nd_order (A, i, j) result (Aij)

    real(spamm_kind) :: Aij
    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    integer, intent(in) :: i, j

    Aij = 0

  end function spamm_get_matrix_2nd_order

  !> Allocate a 2nd order SpAMM matrix.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !!
  !! @return The new matrix.
  function spamm_allocate_matrix_2nd_order (M, N) result (A)

    integer, intent(in) :: M, N
    type(spamm_matrix_2nd_order), pointer :: A

    integer :: K, matrix_N, number_tiles

    allocate(A)
    A%M = M
    A%N = N

    matrix_N = max(M, N)

    K = CEILING(LOG10(DBLE(matrix_N))/LOG10(2D0))

    ! Double check padded size.
    IF(2**K < matrix_N) THEN
      K = K+1
    ENDIF

    IF(K > 0 .AND. 2**(K-1) > matrix_N) THEN
      K = K-1
    ENDIF

    ! Pad matrix to right size.
    A%N_padded = 2**K

    ! Depth starts from 0:
    number_tiles = CEILING(DBLE(A%N_padded)/SPAMM_BLOCK_SIZE)
    A%depth = FLOOR(LOG(DBLE(number_tiles))/LOG(2D0))

    write(*, *) "allocated ", M, "x", N, " matrix"
    write(*, *) "N_padded = ", A%N_padded
    write(*, *) "depth = ", A%depth

  end function spamm_allocate_matrix_2nd_order

end module spamm_management
