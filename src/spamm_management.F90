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
!! @author Nicolas Bock nicolasbock@freeon.org
module spamm_management

  use spamm_types
  use spamm_globals
  use spamm_utilities

#include "spamm_utility_macros.h"

  implicit none

  private

  PUBLIC :: Copy
  PUBLIC :: Delete
  PUBLIC :: New
  PUBLIC :: NewQuNode
  PUBLIC :: SpAMM_Copy_BiTree_2_BiTree_Recur
  PUBLIC :: SpAMM_Copy_QuTree_2_QuTree_Recur
  PUBLIC :: SpAMM_Delete_QuTree_Recur
  public :: get
  public :: set
  public :: reset_counters
  public :: spamm_allocate_matrix_2nd_order
  public :: spamm_allocate_matrix_order_1
  public :: spamm_identity_matrix
  public :: spamm_zero_matrix
  public :: newbinode
  public :: count_nonzero
  public :: absmax

  !> Interface for deep copies of SpAMM objects.
  INTERFACE Copy
    !MODULE PROCEDURE SpAMM_Copy_QuTree_2_QuTree
    MODULE PROCEDURE SpAMM_Copy_QuTree_2_BiTree
    MODULE PROCEDURE SpAMM_Copy_BiTree_2_BiTree
    module procedure spamm_copy_2nd_order_to_2nd_order
    module procedure spamm_copy_2nd_order_to_order_1
    module procedure spamm_copy_order_1_to_order_1
  END INTERFACE

  !> Interface for deletion (deallocation) of SpAMM objects.
  INTERFACE Delete
    MODULE PROCEDURE SpAMM_Delete_QuTree
    MODULE PROCEDURE SpAMM_Delete_BiTree
    module procedure spamm_delete_matrix_order_1
    module procedure spamm_delete_matrix_2nd_order
  END INTERFACE

  !> Interface for creation (allocation) of SpAMM objects.
  INTERFACE New
    MODULE PROCEDURE SpAMM_Allocate_Full_QuTree
    MODULE PROCEDURE SpAMM_Allocate_Full_BiTree
    module procedure spamm_allocate_matrix_order_1
    module procedure spamm_allocate_matrix_2nd_order
  END INTERFACE

  !> Interface for getting single matrix elements of SpAMM objects.
  interface get
    module procedure spamm_get_matrix_order_1
    module procedure spamm_get_matrix_2nd_order
  end interface get

  !> Interface for setting matrix elements.
  interface set
    module procedure spamm_set_matrix_2nd_order
  end interface set

  !> Interface for reset_counters functions.
  interface reset_counters
    module procedure spamm_reset_counters_2nd_order
  end interface reset_counters

  !> Interface for non-zero element counting functions.
  interface count_nonzero
    module procedure count_nonzero_order_1
    module procedure count_nonzero_order_2
  end Interface count_nonzero

  !> Interface for the absmax functions.
  interface absmax
    module procedure absmax_order_2
  end interface absmax

contains

  !> Copy a 1st order matrix: \f$ B \leftarrow A \f$.
  !!
  !! @param A The vector A.
  !! @param B The vector B.
  subroutine spamm_copy_order_1_to_order_1 (A, B)

    type(spamm_matrix_order_1), pointer, intent(in) :: A
    type(spamm_matrix_order_1), pointer, intent(inout) :: B

    if(.not. associated(A)) then
      return
    endif

    if(associated(B)) then
      call delete(B)
    endif

    LOG_DEBUG("copying vector")

    B => null()
    call spamm_allocate_matrix_order_1(A%N, B)
    call spamm_copy_bitree_2_bitree(A%root, B%root)

  end subroutine spamm_copy_order_1_to_order_1

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

    LOG_DEBUG("copying matrix")

    B => null()
    call spamm_allocate_matrix_2nd_order(A%M, A%N, B)
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
      CALL NewBiNode(bC, qA%i_lower, qA%i_upper)
    endif

    !$OMP TASK UNTIED SHARED(qA,bC)
    CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA, bC, j)
    !$OMP END TASK
    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Copy_QuTree_2_BiTree

  !> Copy bitree into another QuTree: @f$ C \leftarrow A @f$.
  !!
  !! @param bA Vector A.
  !! @param bC Vector C.
  SUBROUTINE SpAMM_Copy_BiTree_2_BiTree(bA,bC)

    TYPE(BiTree), POINTER :: bA,bC
    INTEGER               :: Depth

    IF(.NOT.ASSOCIATED(bC))&
      CALL NewBiNode(bC, bA%i_lower, bA%i_upper)
    Depth=0
    !$OMP TASK UNTIED SHARED(bA,bC)
    CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA,bC,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

  END SUBROUTINE SpAMM_Copy_BiTree_2_BiTree

  !> Delete a 1st order matrix.
  !!
  !! @param A The vector.
  subroutine spamm_delete_matrix_order_1 (A)

    type(spamm_matrix_order_1), pointer, intent(inout) :: A

    call spamm_delete_bitree(A%root)
    deallocate(A)

  end subroutine spamm_delete_matrix_order_1

  !> Delete a 2nd order matrix.
  !!
  !! @param A The matrix.
  subroutine spamm_delete_matrix_2nd_order (A)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A

    if(associated(A)) then
      call spamm_delete_qutree(A%root)
      deallocate(A)
    endif

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

  !> Delete a BiTree: @f$ A \leftarrow 0 @f$.
  !!
  !! @param bA The bitree.
  SUBROUTINE SpAMM_Allocate_Full_BiTree(bA, i_lower, i_upper)

    TYPE(BiTree), POINTER :: bA
    integer, intent(in) :: i_lower, i_upper

    CALL NewBiNode(bA, i_lower, i_upper)
    CALL SpAMM_Allocate_Full_BiTree_Recur(bA)

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

    LOG_DEBUG("q: "//to_string(qA))

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

  !> Recursive copy of bitree to bitree
  !!
  !! @param bA Vector A.
  !! @param bC Vector C.
  !! @param depth The tree depth.
  RECURSIVE SUBROUTINE SpAMM_Copy_BiTree_2_BiTree_Recur(bA,bC,Depth)

    TYPE(BiTree), POINTER  :: bA,bC
    INTEGER                :: Depth

    if(.not. associated(bA)) then
      LOG_DEBUG("bA not associated")
      return
    endif

    if(.not. associated(bC)) then
      call NewBiNode(bC, bA%i_lower, bA%i_upper)
    endif

    bC%Norm=bA%Norm
    if(bA%i_upper-bA%i_lower+1 == SPAMM_BLOCK_SIZE .and. allocated(bA%vect)) then
      if(.not. allocated(bC%Vect)) then
        allocate(bC%Vect(SPAMM_BLOCK_SIZE))
      endif
      LOG_DEBUG("copying vect")
      bC%Vect = bA%Vect
    ELSE
      IF(ASSOCIATED(bA%sect1))THEN
        !$OMP TASK UNTIED SHARED(bA,bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA%sect1,bC%sect1,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(bC%sect1))THEN
        !$OMP TASK UNTIED SHARED(bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_BiTree_Recur(bC%sect1)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(bC%sect1)
      ENDIF
      IF(ASSOCIATED(bA%sect2))THEN
        !$OMP TASK UNTIED SHARED(bA,bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA%sect2,bC%sect2,Depth+1)
        !$OMP END TASK
      ELSEIF(ASSOCIATED(bC%sect2))THEN
        !$OMP TASK UNTIED SHARED(bC) &
        !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Delete_BiTree_Recur(bC%sect2)
        !$OMP END TASK
        !$OMP TASKWAIT
        DEALLOCATE(bC%sect2)
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_BiTree_2_BiTree_Recur

  !> Copy a column from a matrix to a vector.
  !!
  !! @param A The matrix.
  !! @param j The column index.
  !! @param V The vector.
  subroutine spamm_copy_2nd_order_to_order_1 (A, j, V)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    integer, intent(in) :: j
    type(spamm_matrix_order_1), pointer, intent(inout) :: V

    if(.not. associated(A)) then
      return
    endif

    if(associated(V)) then
      if(V%N /= A%M) then
        LOG_FATAL("dimension mismatch")
        error stop
      endif
    else
      call spamm_allocate_matrix_order_1(A%M, V)
    endif

    LOG_DEBUG("copying column "//to_string(j)//" from A to V")

    call spamm_copy_qutree_2_bitree_recur(A%root, V%root, j)

  end subroutine spamm_copy_2nd_order_to_order_1

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

    if(.not. associated(qA)) then
      LOG_DEBUG("A not associated")
      return
    endif

    IF(.NOT.ASSOCIATED(bC)) then
      CALL NewBiNode(bC, qA%i_lower, qA%i_upper)
    endif

    LOG_DEBUG("A: "//to_string(qA))
    LOG_DEBUG("V: "//to_string(BC))

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      if(allocated(qA%blok)) then
        IF(.NOT. ALLOCATED(bC%Vect)) then
          ALLOCATE(bC%Vect(SPAMM_BLOCK_SIZE))
        endif
        bC%Vect = qA%Blok(:, j-qA%j_lower+1)
      else
        bC%Vect = SpAMM_Zero
      endif
    ELSE
      half = (qA%j_upper-qA%j_lower+1)/2-1
      LOG_DEBUG("qA%j_lower+half = "//to_string(qA%j_lower+half))
      if(j <= qA%j_lower+half) then
        LOG_DEBUG("descending upper half")
        IF(ASSOCIATED(qA%Quad11))THEN
          !$OMP TASK UNTIED SHARED(qA,bC)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad11, bC%sect1, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%sect1))THEN
          !$OMP TASK UNTIED SHARED(bC)
          CALL SpAMM_Delete_BiTree_Recur(bC%sect1)
          !$OMP END TASK
          DEALLOCATE(bC%sect1)
        ENDIF

        IF(ASSOCIATED(qA%Quad21))THEN
          !$OMP TASK UNTIED SHARED(qA,bC)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad21, bC%sect2, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%sect2))THEN
          !$OMP TASK UNTIED SHARED(bC)
          CALL SpAMM_Delete_BiTree_Recur(bC%sect2)
          !$OMP END TASK
          DEALLOCATE(bC%sect2)
        ENDIF
      ELSE
        LOG_DEBUG("descending lower half")
        IF(ASSOCIATED(qA%Quad12))THEN
          !$OMP TASK UNTIED SHARED(qA,bC)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad12, bC%sect1, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%sect1))THEN
          !$OMP TASK UNTIED SHARED(bC)
          CALL SpAMM_Delete_BiTree_Recur(bC%sect1)
          !$OMP END TASK
          DEALLOCATE(bC%sect1)
        ENDIF

        IF(ASSOCIATED(qA%Quad22))THEN
          !$OMP TASK UNTIED SHARED(qA,bC)
          CALL SpAMM_Copy_QuTree_2_BiTree_Recur(qA%Quad22, bC%sect2, j)
          !$OMP END TASK
        ELSEIF(ASSOCIATED(bC%sect2))THEN
          !$OMP TASK UNTIED SHARED(bC)
          CALL SpAMM_Delete_BiTree_Recur(bC%sect2)
          !$OMP END TASK
          DEALLOCATE(bC%sect2)
        ENDIF
      ENDIF
      !$OMP TASKWAIT
      LOG_DEBUG("going back up")
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
    IF(ASSOCIATED(bA%sect1))THEN
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_BiTree_Recur(bA%sect1)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(bA%sect1)
      !$OMP END CRITICAL
    ENDIF
    IF(ASSOCIATED(bA%sect2))THEN
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_BiTree_Recur(bA%sect2)
      !$OMP END TASK
      !$OMP TASKWAIT
      !$OMP CRITICAL
      DEALLOCATE(bA%sect2)
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

  !> Allocate new bitree node.
  !!
  !! @param bA The node pointer.
  !! @param i_lower The lower index.
  !! @param i_upper The upper index.
  SUBROUTINE NewBiNode(bA, i_lower, i_upper)

    TYPE(BiTree), POINTER :: bA
    integer, intent(in) :: i_lower, i_upper

    if(associated(bA)) then
      call delete(bA)
    endif

    allocate(bA)

    bA%i_lower = i_lower
    bA%i_upper = i_upper

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

    LOG_DEBUG("q: "//to_string(i_lower)//" "//to_string(i_upper))
    LOG_DEBUG("   "//to_string(j_lower)//" "//to_string(j_upper))

    if(rows /= columns) then
      LOG_FATAL("non-square submatrix")
      error stop
    endif

    if(rows < SPAMM_BLOCK_SIZE .or. columns < SPAMM_BLOCK_SIZE) then
      LOG_FATAL("[ozC7x7z3HIgTIa0Q] logic error")
      LOG_FATAL("rows = "//to_string(rows))
      LOG_FATAL("columns = "//to_string(columns))
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

  RECURSIVE SUBROUTINE SpAMM_Allocate_Full_BiTree_Recur(bA)

    TYPE(BiTree), POINTER :: bA
    integer :: half

    IF(bA%i_upper-bA%i_lower+1 == SPAMM_BLOCK_SIZE)THEN
      ALLOCATE(bA%Vect(SPAMM_BLOCK_SIZE))
      bA%Vect=SpAMM_Zero
    ELSE
      half = (bA%i_upper-bA%i_lower+1)/2-1
      call newbinode(bA%sect1, bA%i_lower, bA%i_lower+half)
      call newbinode(bA%sect2, bA%i_lower+half+1, bA%i_upper)
      CALL SpAMM_Allocate_Full_BiTree_Recur(bA%sect1)
      CALL SpAMM_Allocate_Full_BiTree_Recur(bA%sect2)
    ENDIF

  END SUBROUTINE SpAMM_Allocate_Full_BiTree_Recur

  !> Construct the identy matrix of size M x N.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !!
  !! @return The matrix.
  function spamm_identity_matrix (M, N) result (A)

    integer, intent(in) :: M, N
    type(spamm_matrix_2nd_order), pointer :: A
    integer :: i

    if(M /= N) then
      LOG_FATAL("M /= N")
      error stop
    endif

    LOG_DEBUG("creating "//to_string(M)//"x"//to_string(N)//" identity matrix")

    A => null()
    call spamm_allocate_matrix_2nd_order(M, N, A)

    do i = 1, N
      call spamm_set_matrix_2nd_order(A, i, i, 1.0_spamm_kind)
    enddo

  end function spamm_identity_matrix

  !> Construct a zero matrix of size M x N.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !!
  !! @return The matrix.
  function spamm_zero_matrix (M, N) result (A)

    type(spamm_matrix_2nd_order), pointer :: A
    integer, intent(in) :: M, N

    A => null()
    call spamm_allocate_matrix_2nd_order(M, N, A)
    call spamm_allocate_full_qutree_recur(A%root, 1, 1, A%N_padded, A%N_padded)

  end function spamm_zero_matrix

  !> Get a vector element from a 1st order SpAMM matrix.
  !!
  !! @param A The matrix.
  !! @param i The row index.
  !!
  !! @return The matrix element.
  function spamm_get_matrix_order_1 (V, i) result (Vi)

    real(spamm_kind) :: Vi
    type(spamm_matrix_order_1), pointer, intent(in) :: V
    integer, intent(in) :: i

    LOG_DEBUG("getting element "//to_string(i))
    Vi = spamm_get_bitree(V%root, i)

  end function spamm_get_matrix_order_1

  !> Get a vector element from a spamm_types::bitree.
  !!
  !! @param qV A pointer to a bitree.
  !! @param i The row index.
  !!
  !! @return The matrix element.
  recursive function spamm_get_bitree (qV, i) result (Vi)

    real(spamm_kind) :: Vi
    type(bitree), pointer, intent(in) :: qV
    integer, intent(in) :: i

    Vi = 0

    LOG_DEBUG("q: "//to_string(qV%i_lower)//", "//to_string(qV%i_upper))

    if(.not. associated(qV)) return

    if(i > qV%i_upper) then
      LOG_FATAL("logic error, i ("//to_string(i)//") above upper bound ("//to_string(qV%i_upper)//")")
      error stop
    endif

    if(i < qV%i_lower) then
      LOG_FATAL("logic error, i ("//to_string(i)//") below lower bound ("//to_string(qV%i_lower)//")")
      error stop
    endif

    if(qV%i_upper-qV%i_lower+1 == SPAMM_BLOCK_SIZE) then
      if(allocated(qV%vect)) then
        Vi = qV%vect(i-qV%i_lower+1)
        LOG_DEBUG("found matrix element: Vi = "//to_string(Vi))
      endif
    else
      if(associated(qV%sect1)) then
        if(i >= qV%sect1%i_lower .and. i <= qV%sect1%i_upper) then
          Vi = spamm_get_bitree(qV%sect1, i)
          return
        endif
      endif

      if(associated(qV%sect2)) then
        if(i >= qV%sect2%i_lower .and. i <= qV%sect2%i_upper) then
          Vi = spamm_get_bitree(qV%sect2, i)
          return
        endif
      endif
    endif

  end function spamm_get_bitree

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

    if(associated(A)) then
      Aij = spamm_get_qutree(A%root, i, j)
    else
      Aij = 0
    endif

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

    if(i > qA%i_upper) then
      LOG_FATAL("logic error, i ("//to_string(i)//") above upper bound ("//to_string(qA%i_upper)//")")
      error stop
    endif

    if(j > qA%j_upper) then
      LOG_FATAL("logic error, j ("//to_string(j)//") above upper bound ("//to_string(qA%j_upper)//")")
      error stop
    endif

    if(i < qA%i_lower) then
      LOG_FATAL("logic error, i ("//to_string(i)//") below lower bound ("//to_string(qA%i_lower)//")")
      error stop
    endif

    if(j < qA%j_lower) then
      LOG_FATAL("logic error, j ("//to_string(j)//") below lower bound ("//to_string(qA%j_lower)//")")
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

  !> Set a matrix element from a 2nd order SpAMM matrix.
  !!
  !! @param A The matrix.
  !! @param i The row index.
  !! @param j The column index.
  !! @param Aij The matrix element.
  subroutine spamm_set_matrix_2nd_order (A, i, j, Aij)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    integer, intent(in) :: i, j
    real(spamm_kind), intent(in) :: Aij

    LOG_DEBUG("setting matrix element A("//to_string(i)//","//to_string(j)//")")

    if(.not. associated(A)) then
      LOG_FATAL("A is not associated")
      error stop
    endif

    if(.not. associated(A%root)) then
      allocate(A%root)
      A%root%i_lower = 1
      A%root%i_upper = A%N_padded
      A%root%j_lower = 1
      A%root%j_upper = A%N_padded
    endif

    call spamm_set_qutree(A%root, i, j, Aij)

    A%norm = A%root%norm

    LOG_DEBUG("done setting")

  end subroutine spamm_set_matrix_2nd_order

  !> Set a matrix element from a spamm_types::qutree.
  !!
  !! @param qA A pointer to a qutree.
  !! @param i The row index.
  !! @param j The column index.
  !! @param Aij The matrix element.
  recursive subroutine spamm_set_qutree (qA, i, j, Aij)

    type(qutree), pointer, intent(in) :: qA
    integer, intent(in) :: i, j
    real(spamm_kind), intent(in) :: Aij
    integer :: half

    LOG_DEBUG("qutree setting matrix element A("//to_string(i)//","//to_string(j)//")")
    LOG_DEBUG("  "//to_string(qA))

    if(i > qA%i_upper .or. j > qA%j_upper) then
      LOG_FATAL("logic error, i or j above upper bound")
      LOG_FATAL("i = "//to_string(i)//", j = "//to_string(j))
      LOG_FATAL("i_upper = "//to_string(qA%i_upper)//", j_upper = "//to_string(qA%j_upper))
      error stop
    endif

    if(i < qA%i_lower .or. j < qA%j_lower) then
      LOG_FATAL("logic error, i or j below lower bound")
      LOG_FATAL("i = "//to_string(i)//", j = "//to_string(j))
      LOG_FATAL("i_lower = "//to_string(qA%i_lower)//", j_lower = "//to_string(qA%j_lower))
      error stop
    endif

    if(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE .and. qA%j_upper-qA%j_lower+1 == SPAMM_BLOCK_SIZE) then
      if(.not. allocated(qA%blok)) then
        LOG_DEBUG("allocating new blok")
        allocate(qA%blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
        qA%blok = 0
      endif
      LOG_DEBUG("setting blok("//to_string(i-qA%i_lower+1)//","//to_string(j-qA%j_lower+1)//")")
      qA%blok(i-qA%i_lower+1, j-qA%j_lower+1) = Aij
      qA%norm = sqrt(sum(qA%blok**2))
    else
      half = (qA%i_upper-qA%i_lower+1)/2-1

      if(i <= qA%i_lower+half .and. j <= qA%j_lower+half) then
        if(.not. associated(qA%quad11)) then
          LOG_DEBUG("allocating quad11")
          allocate(qA%quad11)
          qA%quad11%i_lower = qA%i_lower
          qA%quad11%i_upper = qA%i_lower+half
          qA%quad11%j_lower = qA%j_lower
          qA%quad11%j_upper = qA%j_lower+half
        endif
        LOG_DEBUG("going down quad11")
        call spamm_set_qutree(qA%quad11, i, j, Aij)
      else if(i <= qA%i_lower+half .and. j >= qA%j_lower+half+1) then
        if(.not. associated(qA%quad12)) then
          LOG_DEBUG("allocating quad12")
          allocate(qA%quad12)
          qA%quad12%i_lower = qA%i_lower
          qA%quad12%i_upper = qA%i_lower+half
          qA%quad12%j_lower = qA%j_lower+half+1
          qA%quad12%j_upper = qA%j_upper
        endif
        LOG_DEBUG("going down quad12")
        call spamm_set_qutree(qA%quad12, i, j, Aij)
      else if(i >= qA%i_lower+half+1 .and. j <= qA%j_lower+half) then
        if(.not. associated(qA%quad21)) then
          LOG_DEBUG("allocating quad21")
          allocate(qA%quad21)
          qA%quad21%i_lower = qA%i_lower+half+1
          qA%quad21%i_upper = qA%i_upper
          qA%quad21%j_lower = qA%j_lower
          qA%quad21%j_upper = qA%j_lower+half
        endif
        LOG_DEBUG("going down quad21")
        call spamm_set_qutree(qA%quad21, i, j, Aij)
      else if(i >= qA%i_lower+half+1 .and. j >= qA%j_lower+half+1) then
        if(.not. associated(qA%quad22)) then
          LOG_DEBUG("allocating quad22")
          allocate(qA%quad22)
          qA%quad22%i_lower = qA%i_lower+half+1
          qA%quad22%i_upper = qA%i_upper
          qA%quad22%j_lower = qA%j_lower+half+1
          qA%quad22%j_upper = qA%j_upper
        endif
        LOG_DEBUG("going down quad22")
        call spamm_set_qutree(qA%quad22, i, j, Aij)
      endif

      qA%norm = 0

      if(associated(qA%quad11)) then
        qA%norm = qA%norm + qA%quad11%norm**2
      endif

      if(associated(qA%quad12)) then
        qA%norm = qA%norm + qA%quad12%norm**2
      endif

      if(associated(qA%quad21)) then
        qA%norm = qA%norm + qA%quad21%norm**2
      endif

      if(associated(qA%quad22)) then
        qA%norm = qA%norm + qA%quad22%norm**2
      endif

      qA%norm = sqrt(qA%norm)
    endif

    LOG_DEBUG("going back up")

  end subroutine spamm_set_qutree

  !> Allocate a 1st order SpAMM matrix.
  !!
  !! @param N The number of columns.
  !! @param V The new vector.
  subroutine spamm_allocate_matrix_order_1 (N, V)

    integer, intent(in) :: N
    type(spamm_matrix_order_1), pointer, intent(out) :: V

    V => null()
    allocate(V)

    V%N = N
    V%depth = 0
    V%N_padded = SPAMM_BLOCK_SIZE

    ! This should be pretty efficient for reasonable matrix sizes and is
    ! presumably faster than some logarithm calculation since it only involves
    ! an increment and a bit shift.
    do while(V%N_padded < N)
      LOG_DEBUG("depth = "//to_string(V%depth)//", N_padded = "//to_string(V%N_padded))
      V%depth = V%depth+1
      V%N_padded = 2*V%N_padded
    enddo

    LOG_INFO("allocated "//to_string(N)//" vector")
    LOG_INFO("  BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE))
    LOG_INFO("  N_padded   = "//to_string(V%N_padded))
    LOG_INFO("  depth      = "//to_string(V%depth))

  end subroutine spamm_allocate_matrix_order_1

  !> Allocate a 2nd order SpAMM matrix.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !! @param A The new matrix.
  subroutine spamm_allocate_matrix_2nd_order (M, N, A)

    integer, intent(in) :: M, N
    type(spamm_matrix_2nd_order), pointer, intent(out) :: A

    A => null()
    allocate(A)

    A%M = M
    A%N = N
    A%depth = 0
    A%N_padded = SPAMM_BLOCK_SIZE

    ! This should be pretty efficient for reasonable matrix sizes and is
    ! presumably faster than some logarithm calculation since it only involves
    ! an increment and a bit shift.
    do while(A%N_padded < max(M, N))
      LOG_DEBUG("depth = "//to_string(A%depth)//", N_padded = "//to_string(A%N_padded))
      A%depth = A%depth+1
      A%N_padded = 2*A%N_padded
    enddo

    LOG_DEBUG("allocated "//to_string(M)//"x"//to_string(N)//" matrix")
    LOG_DEBUG("  BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE))
    LOG_DEBUG("  N_padded   = "//to_string(A%N_padded))
    LOG_DEBUG("  depth      = "//to_string(A%depth))

  end subroutine spamm_allocate_matrix_2nd_order

  !> Reset all counters in a qutree recursively.
  !!
  !! @param qA The qutree node.
  recursive subroutine spamm_reset_counters_qutree (qA)

    type(qutree), pointer, intent(inout) :: qA

    if(.not. associated(qA)) return

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

    LOG_DEBUG("resetting counters")

    A%number_operations = 0
    A%number_nonzeros = 0
    call spamm_reset_counters_qutree(A%root)

  end subroutine spamm_reset_counters_2nd_order

  !> Count the number of non-zero elements in a dense matrix.
  !!
  !! @param A The dense matrix.
  !!
  !! @return The number of non-zero elements.
  function count_nonzero_order_2 (A) result(number_nonzeros)

    real(spamm_kind), intent(in) :: A(:, :)
    real(spamm_kind) :: number_nonzeros

    integer :: i, j

    number_nonzeros = sum(reshape( &
      [ ((1, i = 1, size(A, 1)), j = 1, size(A, 2)) ], &
      [ size(A, 1), size(A, 2) ]), &
      reshape( &
      [ ((A(i, j) /= 0.0, i = 1, size(A, 1)), j = 1, size(A, 2)) ], &
      [ size(A, 1), size(A, 2) ]))

  end function count_nonzero_order_2

  !> Count the number of non-zero elements in a dense vector.
  !!
  !! @param A The dense vector.
  !!
  !! @return The number of non-zero elements.
  function count_nonzero_order_1 (A) result(number_nonzeros)

    real(spamm_kind), intent(in) :: A(:)
    real(spamm_kind) :: number_nonzeros

    integer :: i

    number_nonzeros = sum(reshape( &
      [ (1, i = 1, size(A, 1)) ], [ size(A, 1) ]), &
      reshape((/ (A(i) /= 0.0, i = 1, size(A, 1)) /), (/ size(A, 1) /)))

  end function count_nonzero_order_1

  !> Find the maximum matrix element, \f$ \max | A_{ij} | \f$.
  !!
  !! @param A The matrix node.
  !!
  !! @return The maximum matrix element.
  recursive function absmax_qutree (A) result(absmax)

    type(qutree), pointer, intent(in) :: A
    real(kind(0d0)) :: absmax

    absmax = 0

    if(.not. associated(A)) then
      return
    endif

    if(A%i_upper-A%i_lower+1 == SPAMM_BLOCK_SIZE) then
      absmax = maxval(abs(A%blok))
    else
      if(associated(A%quad11)) then
        absmax = max(absmax, absmax_qutree(A%quad11))
      endif
      if(associated(A%quad12)) then
        absmax = max(absmax, absmax_qutree(A%quad12))
      endif
      if(associated(A%quad21)) then
        absmax = max(absmax, absmax_qutree(A%quad21))
      endif
      if(associated(A%quad22)) then
        absmax = max(absmax, absmax_qutree(A%quad22))
      endif
    endif

  end function absmax_qutree

  !> Find the maximum matrix element, \f$ \max | A_{ij} | \f$.
  !!
  !! @param A The matrix.
  !!
  !! @return The maximum matrix element.
  function absmax_order_2 (A) result(absmax)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    real(spamm_kind) :: absmax

    absmax = 0

    if(.not. associated(A)) then
      return
    endif

    absmax = absmax_qutree(A%root)

  end function absmax_order_2

end module spamm_management
