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
  public :: reset_counters

  !> Interface for deep copies of SpAMM objects.
  INTERFACE Copy
    !MODULE PROCEDURE SpAMM_Copy_QuTree_2_QuTree
    MODULE PROCEDURE SpAMM_Copy_QuTree_2_BiTree
    MODULE PROCEDURE SpAMM_Copy_BiTree_2_BiTree
    module procedure spamm_copy_2nd_order_to_2nd_order
  END INTERFACE

  !> Interface for deletion (deallocation) of SpAMM objects.
  INTERFACE Delete
    MODULE PROCEDURE SpAMM_Delete_QuTree
    MODULE PROCEDURE SpAMM_Delete_BiTree
    module procedure spamm_delete_matrix_2nd_order
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

  interface reset_counters
    module procedure spamm_reset_counters_2nd_order
  end interface reset_counters

CONTAINS

  !> Copy a 2nd order matrix: \f$ B \leftarrow A \f$.
  !!
  !! @param A The matrix A.
  !! @param B The matrix B.
  subroutine spamm_copy_2nd_order_to_2nd_order (A, B)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    type(spamm_matrix_2nd_order), pointer, intent(inout) :: B

    if(.not. associated(A)) then
      return
    endif

    if(associated(B)) then
      call delete(B)
    endif

    B => spamm_allocate_matrix_2nd_order(A%M, A%N)
    call spamm_copy_qutree_2_qutree_recur(A%root, B%root)

  end subroutine spamm_copy_2nd_order_to_2nd_order

  !> Copy a coloumn of A into vector C: \f$ C \leftarrow A(:, j) \f$.
  !!
  !! @param qA Pointer to quadtree.
  !! @param j The column index.
  !! @param bC Pointer to bitree.
  SUBROUTINE SpAMM_Copy_QuTree_2_BiTree (qA, j, bC)

    TYPE(QuTree), POINTER, intent(in) :: qA
    TYPE(BiTree), POINTER, intent(inout) :: bC
    INTEGER, intent(in) :: j

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

  !> Delete a 2nd order matrix.
  !!
  !! @param A The matrix.
  subroutine spamm_delete_matrix_2nd_order (A)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A

    call spamm_delete_qutree(A%root)
    deallocate(A)

  end subroutine spamm_delete_matrix_2nd_order

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

  !> Recursive copy of quadtree: \f$ C \leftarrow A \f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qC Pointer to matrix C.
  RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur (qA, qC)

    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC

    IF(.NOT.ASSOCIATED(qA)) RETURN

    IF(.NOT.ASSOCIATED(qC)) THEN
      CALL NewQuNode(qC, qA%i_lower, qA%j_lower, qA%i_upper, qA%j_upper)
    ENDIF

#if DEBUG >= 2
    write(*, *) "q: ", qA%i_lower, qA%i_upper, qA%j_lower, qA%j_upper
#endif

    qC%Norm = qA%Norm

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      IF(.NOT. allocated(qC%Blok)) THEN
        ALLOCATE(qC%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
      ENDIF
      qC%Blok = qA%Blok
    ELSE
      IF(ASSOCIATED(qA%Quad11)) THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad11, qC%Quad11)
        !$OMP END TASK
      ENDIF
      IF(ASSOCIATED(qA%Quad12)) THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad12, qC%Quad12)
        !$OMP END TASK
      ENDIF
      IF(ASSOCIATED(qA%Quad21)) THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad21, qC%Quad21)
        !$OMP END TASK
      ENDIF
      IF(ASSOCIATED(qA%Quad22)) THEN
        !$OMP TASK UNTIED SHARED(qA,qC)
        CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad22, qC%Quad22)
        !$OMP END TASK
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

    INTEGER               :: half
    INTEGER,DIMENSION(2)  :: Cols,Col_00_10,Col_01_11

    integer :: depth

    depth = 0

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
    INTEGER               :: Status, Depth

    depth = 0

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

    qA%i_lower = i_lower
    qA%i_upper = i_upper
    qA%j_lower = j_lower
    qA%j_upper = j_upper

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
      IF(ASSOCIATED(bA))STOP '[PdMh8CTirrvi44hv] LOGIC ERROR IN NewBiNode'
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

#if DEBUG >= 2
    write(*, *) "q: ", i_lower, i_upper, j_lower, j_upper
#endif

    if(rows /= columns) then
      write(*, *) "non-square submatrix"
      error stop
    endif

    if(rows < SPAMM_BLOCK_SIZE .or. columns < SPAMM_BLOCK_SIZE) then
      write(*, *) "[ozC7x7z3HIgTIa0Q] logic error"
      write(*, *) "rows = ", rows
      write(*, *) "columns = ", columns
      error stop
    endif

    ! Allocate new node.
    CALL NewQuNode(qA, i_lower, j_lower, i_upper, j_upper)

    IF(rows == SPAMM_BLOCK_SIZE) then
      ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
      qA%Blok = SpAMM_Zero
      RETURN
    ELSE
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad11, i_lower, j_lower, i_lower+rows/2-1, j_lower+columns/2-1)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad12, i_lower, j_lower+columns/2, i_lower+rows/2-1, j_upper)
      CALL SpAMM_Allocate_Full_QuTree_Recur(qA%Quad21, i_lower+rows/2, j_lower, i_upper, j_lower+columns/2-1)
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
    call spamm_allocate_full_qutree_recur(A%root, 1, 1, A%N_padded, A%N_padded)

  end function spamm_zero_matrix

  !> Get a matrix element from a 2nd order SpAMM matrix.
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

    Aij = spamm_get_qutree(A%root, i, j)

  end function spamm_get_matrix_2nd_order

  !> Get a matrix element from a spamm_types::qutree.
  !!
  !! @param qA A pointer to a qutree.
  !! @param i The row index.
  !! @param j The column index.
  !!
  !! @return The matrix element.
  recursive function spamm_get_qutree (qA, i, j) result (Aij)

    real(spamm_kind) :: Aij
    type(qutree), pointer, intent(in) :: qA
    integer, intent(in) :: i, j

    Aij = 0

    if(.not. associated(qA)) return

    if(i > qA%i_upper .or. j > qA%j_upper) then
      write(*, *) "[F8xYAsM46GYfJP2j] logic error, i or j above upper bound"
      error stop
    endif

    if(i < qA%i_lower .or. j < qA%j_lower) then
      write(*, *) "[3lJNYprqCQWU3ACZ] logic error, i or j below lower bound"
      error stop
    endif

    if(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE .and. qA%j_upper-qA%j_lower+1 == SPAMM_BLOCK_SIZE) then
      if(allocated(qA%blok)) then
        Aij = qA%blok(i-qA%i_lower+1, j-qA%j_lower+1)
      endif
    else
      if(associated(qA%quad11)) then
        if(i >= qA%quad11%i_lower .and. i <= qA%quad11%i_upper .and. j >= qA%quad11%j_lower .and. j <= qA%quad11%j_upper) then
          Aij = spamm_get_qutree(qA%quad11, i, j)
          return
        endif
      endif

      if(associated(qA%quad12)) then
        if(i >= qA%quad12%i_lower .and. i <= qA%quad12%i_upper .and. j >= qA%quad12%j_lower .and. j <= qA%quad12%j_upper) then
          Aij = spamm_get_qutree(qA%quad12, i, j)
          return
        endif
      endif

      if(associated(qA%quad21)) then
        if(i >= qA%quad21%i_lower .and. i <= qA%quad21%i_upper .and. j >= qA%quad21%j_lower .and. j <= qA%quad21%j_upper) then
          Aij = spamm_get_qutree(qA%quad21, i, j)
          return
        endif
      endif

      if(associated(qA%quad22)) then
        if(i >= qA%quad22%i_lower .and. i <= qA%quad22%i_upper .and. j >= qA%quad22%j_lower .and. j <= qA%quad22%j_upper) then
          Aij = spamm_get_qutree(qA%quad22, i, j)
          return
        endif
      endif
    endif

  end function spamm_get_qutree

  !> Allocate a 2nd order SpAMM matrix.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !!
  !! @return The new matrix.
  function spamm_allocate_matrix_2nd_order (M, N) result (A)

    integer, intent(in) :: M, N
    type(spamm_matrix_2nd_order), pointer :: A

    allocate(A)
    A%M = M
    A%N = N

    A%depth = 0
    A%N_padded = SPAMM_BLOCK_SIZE

    ! This should be pretty efficient for reasonable matrix sizes and is presumably faster than some logarithm calculation since it
    ! only involves an increment and a bit shift.
    do while(A%N_padded < max(M, N))
#if DEBUG >= 1
      write(*, *) "depth = ", A%depth, ", N_padded = ", A%N_padded
#endif
      A%depth = A%depth+1
      A%N_padded = 2*A%N_padded
    enddo

#if DEBUG >= 1
    write(*, *) "[allocate 2nd-order] allocated    ", M, "x", N, " matrix"
    write(*, *) "[allocate 2nd-order] BLOCK_SIZE = ", SPAMM_BLOCK_SIZE
    write(*, *) "[allocate 2nd-order] N_padded   = ", A%N_padded
    write(*, *) "[allocate 2nd-order] depth      = ", A%depth
#endif

  end function spamm_allocate_matrix_2nd_order

  !> Reset all counters in a qutree recursively.
  !!
  !! @param qA The qutree node.
  recursive subroutine spamm_reset_counters_qutree (qA)

    type(qutree), pointer, intent(inout) :: qA

    qA%number_operations = 0
    qA%number_nonzeros = 0

    if(associated(qA%quad11)) then
      call spamm_reset_counters_qutree(qA%quad11)
    endif

    if(associated(qA%quad12)) then
      call spamm_reset_counters_qutree(qA%quad12)
    endif

    if(associated(qA%quad21)) then
      call spamm_reset_counters_qutree(qA%quad21)
    endif

    if(associated(qA%quad22)) then
      call spamm_reset_counters_qutree(qA%quad22)
    endif

  end subroutine spamm_reset_counters_qutree

  !> Reset all counters in a 2nd order SpAMM matrix.
  !!
  !! @param A The matrix.
  subroutine spamm_reset_counters_2nd_order (A)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A

    A%number_operations = 0
    A%number_nonzeros = 0
    call spamm_reset_counters_qutree(A%root)

  end subroutine spamm_reset_counters_2nd_order

end module spamm_management
