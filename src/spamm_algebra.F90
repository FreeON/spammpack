!> Defines operations on SpAMM trees.
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
MODULE SpAMM_ALGEBRA

  use spamm_types
  use spamm_globals
  use spamm_management

#ifdef _OPENMP
  use omp_lib
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Multiply
  PUBLIC :: Trace
  PUBLIC :: Add
  PUBLIC :: Filter
  PUBLIC :: Norm
  PUBLIC :: Dot

  !> The constant @f$ \alpha @f$ for adding 2 quadtrees in place.
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_InPlace_Beta
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_RePlace_Beta
  REAL(SpAMM_KIND) :: SpAMM_Add_Identity_2_QuTree_InPlace_Alpha

  !> The SpAMM threshold for the multiply.
  REAL(SpAMM_KIND) :: SpAMM_Threshold_Multiply_QuTree_x_QuTree
  REAL(SpAMM_KIND) :: SpAMM_Threshold_Multiply_QuTree_x_BiTree

  !> Interface for multiplication operations between different SpAMM types.
  !!
  !! The Sparse Approximate Matrix-Multiply (SpAMM):
  !! @f$ C \leftarrow A \times B @f$.
  INTERFACE Multiply
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_QuTree
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_Scalar
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_BiTree
    MODULE PROCEDURE SpAMM_Multiply_BiTree_x_Scalar
    module procedure spamm_multiply_2nd_order_x_2nd_order
    module procedure spamm_multiply_2nd_order_x_scalar
  END INTERFACE

  !> Interface for trace operations.
  INTERFACE Trace
    MODULE PROCEDURE SpAMM_Trace_QuTree
    MODULE PROCEDURE SpAMM_Trace_QuTree_Product
  END INTERFACE

  !> Interface for additions operations between different SpAMM types.
  INTERFACE Add
    MODULE PROCEDURE SpAMM_Add_QuTree_2_QuTree_InPlace
    MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_InPlace
    MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_RePlace
    MODULE PROCEDURE SpAMM_Add_Identity_2_QuTree_InPlace
    module procedure spamm_add_identity_to_matrix_2nd_order
  END INTERFACE

  !> Interface for filter operations (thresholding of small matrix elements).
  INTERFACE Filter
    MODULE PROCEDURE SpAMM_Filter_QuTree
  END INTERFACE

  !> Interface for norm operations.
  INTERFACE Norm
    MODULE PROCEDURE SpAMM_Norm_Reduce_BiTree
    MODULE PROCEDURE SpAMM_Norm_Reduce_QuTree
    module procedure spamm_norm_reduce_matrix_2nd_order
  END INTERFACE

  !> Interface for dot product operations.
  INTERFACE Dot
    MODULE PROCEDURE SpAMM_Dot_Product_BiTree
  END INTERFACE

CONTAINS

  !> Multiplication operation between a quadtree and a quadtree, @f$ C \leftarrow A \times B @f$.
  !!
  !! @param qA Pointer to quadtree A.
  !! @param qB Pointer to quadtree B.
  !! @param qC Pointer to quadtree C.
  !! @param threshold The SpAMM threshold overriding the global value, spamm_types::spamm_product_tolerance.
  SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree(qA, qB, qC, threshold)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA, qB
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
    REAL(SpAMM_KIND), OPTIONAL :: threshold

    real(spamm_kind) :: local_threshold
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()

    !$OMP PARALLEL

    ! The master thread will lead execution of the product. All subsequent tasks
    ! are untied and can be executed by any thread in the thread group.
    !$OMP MASTER

#ifdef _OPENMP
    WRITE(*,*) "Multiply on ", omp_get_num_threads(), " OpenMP threads"
#endif

    IF(PRESENT(threshold))THEN
      local_threshold = threshold
    ELSE
      local_threshold = 0
    ENDIF

    qC%number_operations = 0

    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, SpAMM_Zero)

    !$OMP TASK UNTIED SHARED(qA,qB,qC)
    CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC, qA, qB, local_threshold)
    !$OMP END TASK

    !$OMP END MASTER

    !$OMP END PARALLEL

    qC%norm = norm(qC)
    qC%norm = sqrt(qC%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree

  !> Scalar multiply of 2nd order matrix.
  !!
  !! @param A The matrix.
  !! @param alpha The scalar \f$ \alpha \f$.
  subroutine spamm_multiply_2nd_order_x_scalar (A, alpha)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A
    real(spamm_kind), intent(in) :: alpha

    call spamm_multiply_qutree_x_scalar(A%root, alpha)

  end subroutine spamm_multiply_2nd_order_x_scalar

  !> Scalar multiply: @f$ A \leftarrow a A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param a Scalar a.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar(qA, a)

    TYPE(QuTree), POINTER        :: qA
    REAL(SpAMM_KIND), INTENT(IN) :: a

    INTEGER            :: Depth
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    IF(.NOT.ASSOCIATED(qA))RETURN

    Depth=0
    TInitial = SpAMM_Get_Time()

    !$OMP TASK SHARED(qA)
    CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, a)
    !$OMP END TASK

    !$OMP TASKWAIT

    qA%norm = norm(qA)
    qA%norm = sqrt(qA%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_Scalar",3)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar

  !> Add 2 quadtree matrices, @f$ A \leftarrow \alpha A + \beta B @f$.
  !!
  !! If any of the two factors Alpha or Beta are not supplied then they default
  !! to one.
  !!
  !! @param qA [inout] Pointer to matrix A.
  !! @param qB [in] Pointer to matrix B.
  !! @param Alpha Factor @f$ \alpha @f$.
  !! @param Beta Factor @f$ \beta @f$.
  SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace(qA,qB,Alpha,Beta)

    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    TYPE(QuTree), POINTER, INTENT(IN)    :: qB
    REAL(SpAMM_KIND), OPTIONAL           :: Alpha, Beta
    REAL(SpAMM_KIND)                     :: InPlace_Alpha, InPlace_Beta
    INTEGER                              :: Depth
    REAL(SpAMM_DOUBLE)                   :: TInitial, TTotal

    Depth=0
    IF(PRESENT(Alpha))THEN
      InPlace_Alpha=Alpha
    ELSE
      InPlace_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
      InPlace_Beta=Beta
    ELSE
      InPlace_Beta=SpAMM_One
    ENDIF
    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(qA,qB)
    CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA, qB, &
      InPlace_Alpha, InPlace_Beta, Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

    qA%norm = norm(qA)
    qA%norm = sqrt(qA%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_QuTree_2_QuTree_InPlace",4)

  END SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace

  !> In-place add for 2nd order SpAMM matrix.
  !!
  !! @param alpha The factor \f$ \alpha \f$.
  subroutine spamm_add_identity_to_matrix_2nd_order (A, alpha)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A
    real(spamm_kind), intent(in) :: alpha

    call spamm_add_identity_2_qutree_inplace_recur(A%root, alpha, A%M, A%N)

  end subroutine spamm_add_identity_to_matrix_2nd_order

  !> QuTree In Place Add: \f$ A \leftarrow A + \alpha I \f$.
  !!
  !! @param qA A pointer to a quadtree.
  !! @param alpha The factor \f$ \alpha \f$.
  !! @param M The number of rows of the matrix.
  !! @param N The number of columns of the matrix.
  SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace(qA, alpha, M, N)

    TYPE(QuTree), POINTER, intent(inout) :: qA
    REAL(SpAMM_KIND), intent(in) :: Alpha
    integer, intent(in) :: M, N
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA, alpha, M, N)
    !$OMP END TASK
    !$OMP TASKWAIT

    qA%norm = norm(qA)
    qA%norm = sqrt(qA%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_Identity_2_QuTree_InPlace",5)

  END SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace

  !> Trace for QuTree: @f$ a = \mathrm{Tr} A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !!
  !! @return The trace.
  FUNCTION SpAMM_Trace_QuTree (qA) RESULT(a)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA
    REAL(SpAMM_KIND)                  :: a

    INTEGER            :: Depth
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    Depth=0

    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(qA,a)
    a = SpAMM_Trace_QuTree_Recur(qA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree",6)

  END FUNCTION SpAMM_Trace_QuTree

  !> Trace for quadtree product, @f$ \mathrm{Tr} \left[ A \times B \right] @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qB Pointer to matrix B.
  !!
  !! @return The trace of the matrix produce, @f$ \mathrm{Tr} \left[ A \times B
  !! \right] @f$.
  FUNCTION SpAMM_Trace_QuTree_Product(qA, qB, threshold) RESULT(a)

    TYPE(QuTree), POINTER      :: qA,qB
    REAL(SpAMM_KIND), OPTIONAL :: threshold
    REAL(SpAMM_KIND)           :: a

    INTEGER            :: Depth
    REAL(SpAMM_KIND)   :: multiplyThreshold
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    Depth=0

    IF(PRESENT(threshold)) THEN
      multiplyThreshold = threshold
    ELSE
      multiplyThreshold = SpAMM_PRODUCT_TOLERANCE
    ENDIF

    TInitial = SpAMM_Get_Time()

    !$OMP TASK UNTIED SHARED(qA,a)
    a = SpAMM_Trace_QuTree_Product_Recur(qA, qB, multiplyThreshold, Depth)
    !$OMP END TASK
    !$OMP TASKWAIT

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree_Product",7)

  END FUNCTION SpAMM_Trace_QuTree_Product

  !> Filter for the QuTree: \f$ \tilde{A} = \mathrm{filter} [A, \tau] \f$.
  SUBROUTINE SpAMM_Filter_QuTree(qA,Tau)

    TYPE(QuTree), POINTER  :: qA
    REAL(SpAMM_KIND)       :: Tau
    INTEGER                :: Depth
    REAL(SpAMM_DOUBLE)     :: TInitial, TTotal
    Depth=0
    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Filter_QuTree_Recur(qA,Tau,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Filter_QuTree",8)

  END SUBROUTINE SpAMM_Filter_QuTree

  !> Update max and L2-norm of a matrix.
  !!
  !! @param qA Pointer to matrix.
  !!
  !! @return The norm.
  FUNCTION SpAMM_Norm_Reduce_QuTree(qA) RESULT(Norm)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA
    REAL(SpAMM_KIND) :: Norm
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()
    !$OMP TASK SHARED(Norm,qA)
    Norm = SpAMM_Norm_Reduce_QuTree_Recur(qA)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Norm_Reduce_QuTree",9)

  END FUNCTION SpAMM_Norm_Reduce_QuTree

  !> Add to binary trees: \f$ C \leftarrow \alpha A + \beta B \f$.
  SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace(bC,Alpha,bA,Beta,bB)

    TYPE(BiTree), POINTER, INTENT(IN)    :: bA,bB
    TYPE(BiTree), POINTER, INTENT(INOUT) :: bC
    REAL(SpAMM_KIND), OPTIONAL           :: Alpha,Beta
    INTEGER                              :: Depth
    REAL(SpAMM_DOUBLE)                   :: TInitial, TTotal
    Depth=0
    IF(PRESENT(Alpha))THEN
      SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha=Alpha
    ELSE
      SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
      SpAMM_Add_BiTree_2_BiTree_RePlace_Beta=Beta
    ELSE
      SpAMM_Add_BiTree_2_BiTree_RePlace_Beta=SpAMM_One
    ENDIF
    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(bA,bB)
    CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA,bB,bC,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_BiTree_2_BiTree_RePlace",10)

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace

  !> \f$ A \leftarrow \alpha A + \beta B \f$.
  !!
  !! @param alpha The parameter \f$ \alpha \f$. If not given, \f$ \alpha = 0 \f$.
  !! @param A Vector \f$ A \f$.
  !! @param beta The parameter \f$ \beta \f$. If not given, \f$ \beta = 1 \f$.
  !! @param B Vector \f$ B \f$.
  SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace(bA,Alpha,bB,Beta)

    TYPE(BiTree), POINTER, INTENT(INOUT) :: bA
    TYPE(BiTree), POINTER, INTENT(IN)    :: bB
    REAL(SpAMM_KIND), OPTIONAL           :: Alpha,Beta
    INTEGER                              :: Depth
    REAL(SpAMM_DOUBLE)                   :: TInitial, TTotal

    Depth=0
    IF(PRESENT(Alpha))THEN
      SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha=Alpha
    ELSE
      SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
      SpAMM_Add_BiTree_2_BiTree_InPlace_Beta=Beta
    ELSE
      SpAMM_Add_BiTree_2_BiTree_InPlace_Beta=SpAMM_One
    ENDIF
    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(bA,bB)
    CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA,bB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_BiTree_2_BiTree_InPlace",11)

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace

  !> Inner product: \f$ (A,B) \f$.
  FUNCTION SpAMM_Dot_Product_BiTree(bA,bB) RESULT(dot)

    INTEGER              :: Depth
    TYPE(BiTree),POINTER :: bA,bB
    REAL(SpAMM_KIND)     :: dot
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    TInitial = SpAMM_Get_Time()
    !$OMP TASK SHARED(dot,bA,bB)
    dot=SpAMM_Dot_Product_BiTree_Recur(bA,bB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Dot_Product_BiTree",12)

  END FUNCTION SpAMM_Dot_Product_BiTree

  !> Sparse Approximate Matrix-Vector Multiply (SpAMV): \f$ C \leftarrow A.B \f$.
  SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree(qA,bB,bC,LocalThreshold)

    TYPE(QuTree), POINTER     :: qA
    TYPE(BiTree), POINTER     :: bB,bC
    INTEGER                   :: Depth
    REAL(SpAMM_KIND),OPTIONAL :: LocalThreshold
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    IF(PRESENT(LocalThreshold))THEN
      SpAMM_Threshold_Multiply_QuTree_x_BiTree=LocalThreshold
    ELSE
      SpAMM_Threshold_Multiply_QuTree_x_BiTree=SpAMM_PRODUCT_TOLERANCE
    ENDIF
    Depth=0
    TInitial = SpAMM_Get_Time()
    CALL SpAMM_Multiply_BiTree_x_Scalar(bC,SpAMM_Zero)
    !$OMP TASK UNTIED SHARED(qA,bB,bC)
    CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC,qA,bB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_BiTree",13)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree

  !> Scalar multiply: \f$ A \leftarrow \alpha A \f$.
  RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar(bA,a)

    INTEGER :: Depth
    TYPE(BiTree), POINTER    :: bA
    REAL(SpAMM_KIND)         :: a
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    IF(.NOT.ASSOCIATED(bA))RETURN
    Depth=0
    TInitial = SpAMM_Get_Time()
    !$OMP TASK SHARED(bA)
    CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA,a,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_BiTree_x_Scalar",14)

  END SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar

  !> Norm for BiTrees.
  FUNCTION SpAMM_Norm_Reduce_BiTree(bA) RESULT(Norm)

    INTEGER              :: Depth
    TYPE(BiTree),POINTER :: bA
    REAL(SpAMM_KIND)     :: Norm
    REAL(SpAMM_DOUBLE)   :: TInitial, TTotal

    Depth=0
    TInitial = SpAMM_Get_Time()
    !$OMP TASK SHARED(Norm,bA)
    Norm = SpAMM_Norm_Reduce_BiTree_Recur(bA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Norm_Reduce_BiTree",15)

  END FUNCTION SpAMM_Norm_Reduce_BiTree

  !> Recursive part of the multiplication operation between a quadtree and a
  !! quadtree, @f$ C \leftarrow A \times B @f$.
  !!
  !! @param qA Pointer to quadtree A.
  !! @param qB Pointer to quadtree B.
  !! @param qC Pointer to quadtree C.
  !! @param threshold The SpAMM product tolerance.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur(qC, qA, qB, threshold)

    TYPE(QuTree), POINTER :: qC, qA, qB
    REAL(SpAMM_KIND) :: threshold
    integer :: tier
    INTEGER :: Depth
    REAL(SpAMM_KIND), DIMENSION(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE) :: temp
    INTEGER :: i, j, k, l

    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB)) THEN
      ! Apply the SpAMM condition.
      IF(qA%Norm*qB%Norm < threshold) RETURN
#ifdef _OPENMP
      CALL OMP_SET_LOCK(qC%lock)
#endif
      IF(.NOT.ASSOCIATED(qC))THEN
        ! Allocate new node.
        CALL NewQuNode(qC, qA%i_lower, qB%j_lower, qA%i_upper, qB%j_upper)
      ENDIF
#ifdef _OPENMP
      CALL OMP_UNSET_LOCK(qC%lock)
#endif
      ! At the bottom, calculate the product.
      IF(qC%i_upper-qC%i_lower+1 == SPAMM_BLOCK_SIZE) then
#ifdef _OPENMP
        CALL OMP_SET_LOCK(qC%lock)
#endif
        IF(.NOT.ALLOCATED(qC%Blok))THEN
          ! Allocate new block.
          ALLOCATE(qC%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
          qC%Blok = SpAMM_Zero
        ENDIF
#ifdef _OPENMP
        CALL OMP_UNSET_LOCK(qC%lock)
#endif

#if defined(_OPENMP) && ! defined(BIGLOCK)
        CALL OMP_SET_LOCK(qC%lock)
#endif
        qC%Blok = qC%Blok + MATMUL(qA%Blok, qB%Blok)
        qC%number_operations = qC%number_operations+SPAMM_BLOCK_SIZE**3

#if defined(_OPENMP) && ! defined(BIGLOCK)
        CALL OMP_UNSET_LOCK(qC%lock)
#endif
      ELSE

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad11%Norm*qB%Quad11%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad11, qA%Quad11, qB%Quad11, threshold)
        !$OMP END TASK

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad11%Norm*qB%Quad12%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad12, qA%Quad11, qB%Quad12, threshold)
        !$OMP END TASK

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad21%Norm*qB%Quad11%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad21, qA%Quad21, qB%Quad11, threshold)
        !$OMP END TASK

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad21%Norm*qB%Quad12%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad22, qA%Quad21, qB%Quad12, threshold)
        !$OMP END TASK

#ifdef BIGLOCK
        !$OMP TASKWAIT
#endif

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad12%Norm*qB%Quad21%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad11, qA%Quad12, qB%Quad21, threshold)
        !$OMP END TASK

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad12%Norm*qB%Quad22%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad12, qA%Quad12, qB%Quad22, threshold)
        !$OMP END TASK

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad21%Norm*qB%Quad21%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad21, qA%Quad22, qB%Quad21, threshold)
        !$OMP END TASK

        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(qA%Quad22%Norm*qB%Quad22%Norm > SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad22, qA%Quad22, qB%Quad22, threshold)
        !$OMP END TASK

        !$OMP TASKWAIT

        qC%number_operations = &
          qC%Quad11%number_operations + &
          qC%Quad12%number_operations + &
          qC%Quad21%number_operations + &
          qC%Quad22%number_operations

      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur

  !> Recursive part of scalar multiply with quadtree matrix, @f$ A \leftarrow \alpha A @f$.
  !!
  !! @param qA Pointer to quadtree.
  !! @param alpha The scalar.
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)

    TYPE(QuTree), POINTER :: qA
    REAL(SpAMM_KIND) :: alpha

    IF(.NOT.ASSOCIATED(qA)) RETURN

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      ! At the bottom, multiply the block.
      qA%Norm = qA%Norm*ABS(alpha)
      qA%Blok = alpha*qA%Blok
    ELSE
      !$OMP TASK UNTIED SHARED(qA)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad11, alpha)
      !$OMP END TASK
      !$OMP TASK UNTIED  SHARED(qA)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad12, alpha)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad21, alpha)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad22, alpha)
      !$OMP END TASK
      !$OMP TASKWAIT
      qA%Norm = qA%Norm*ABS(alpha)
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur

  !> Add 2 quadtree matrices, @f$ A \leftarrow \alpha A + \beta B @f$.
  !!
  !! @param qA [inout] Pointer to matrix A.
  !! @param qB [in] Pointer to matrix B.
  !! @param alpha Factor @f$ \alpha @f$.
  !! @param beta Factor @f$ \beta @f$.
  !! @param Depth The current tier.
  RECURSIVE SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA, qB, alpha, beta, Depth)

    TYPE(QuTree),POINTER :: qA,qB
    REAL(SpAMM_KIND)     :: alpha, beta
    INTEGER              :: Depth
    LOGICAL              :: TA, TB

    TA=ASSOCIATED(qA)
    TB=ASSOCIATED(qB)

    IF(TA.AND.TB)THEN
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        qA%Blok=alpha*qA%Blok+beta*qB%Blok
      ELSE
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad11, qB%Quad11, alpha, beta, Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad12, qB%Quad12, alpha, beta, Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad21, qB%Quad21, alpha, beta, Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad22, qB%Quad22, alpha, beta, Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
      ENDIF
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qB, qA)
      !$OMP END TASK
      !$OMP TASKWAIT
    ELSEIF(TA .AND. .NOT.TB) THEN
      ! Multiply A tree with alpha.
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)
    ENDIF

  END SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur

  !> QuTree In Place Add: \f$ A \leftarrow A + \alpha I \f$.
  !!
  !! @param qA A pointer to a quadtree.
  !! @param alpha The factor \f$ \alpha \f$.
  !! @param M The number of rows of the matrix.
  !! @param N The number of columns of the matrix.
  RECURSIVE SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA, alpha, M, N)

    TYPE(QuTree), POINTER, intent(inout) :: qA
    real(spamm_kind), intent(in) :: alpha
    integer, intent(in) :: M, N
    integer :: i, j

    write(*, *) "q:", qA%i_lower, qA%j_lower, qA%i_upper, qA%j_upper

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      if(.not. allocated(qA%blok)) then
        ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
        qA%Blok = SpAMM_Zero
      ENDIF

      write(*, *) "before"
      do i = 1, SPAMM_BLOCK_SIZE
        write(*, "(4f10.3)") (qA%blok(i, j), j = 1, SPAMM_BLOCK_SIZE)
      enddo

      do i = 1, MIN(SPAMM_BLOCK_SIZE, M-qA%i_lower+1, N-qA%j_lower+1)
        qA%Blok(i:i, i:i) = qA%Blok(i:i, i:i)+alpha
      enddo

      write(*, *) "after"
      do i = 1, SPAMM_BLOCK_SIZE
        write(*, "(4f10.3)") (qA%blok(i, j), j = 1, SPAMM_BLOCK_SIZE)
      enddo
    ELSE
      if(associated(qA%quad11)) then
        !$OMP TASK UNTIED SHARED(qA)
        CALL SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA%Quad11, alpha, M, N)
        !$OMP END TASK
      endif

      if(associated(qA%quad22)) then
        !$OMP TASK UNTIED SHARED(qA)
        CALL SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA%Quad22, alpha, M, N)
        !$OMP END TASK
      endif

      !$OMP TASKWAIT
    ENDIF

  END SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace_Recur

  !> Recursive part of trace for QuTree: @f$ a = \mathrm{Tr} A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param Depth The current tier.
  !!
  !! @return The trace.
  RECURSIVE FUNCTION SpAMM_Trace_QuTree_Recur(qA,Depth) RESULT(Trace)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA
    INTEGER, INTENT(IN)               :: Depth
    REAL(SpAMM_KIND)                  :: Trace

    REAL(SpAMM_KIND) :: Trace00, Trace11
    INTEGER          :: I

    Trace = SpAMM_Zero

    IF(Depth == SpAMM_TOTAL_DEPTH) THEN
      Trace = SpAMM_Zero
      IF(.NOT.ASSOCIATED(qA)) RETURN
      DO I = 1, SPAMM_BLOCK_SIZE
        Trace = Trace+qA%Blok(I,I)
      ENDDO
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11).AND. &
        .NOT.ASSOCIATED(qA%Quad22))THEN
      Trace = SpAMM_Zero
    ELSEIF(.NOT.ASSOCIATED(qA%Quad22))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace = SpAMM_Trace_QuTree_Recur(qA%Quad11, Depth+1)
      !$OMP END TASK
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace = SpAMM_Trace_QuTree_Recur(qA%Quad22, Depth+1)
      !$OMP END TASK
    ELSE
      !$OMP TASK UNTIED SHARED(qA,Trace00) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace00 = SpAMM_Trace_QuTree_Recur(qA%Quad11, Depth+1)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,Trace11) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace11 = SpAMM_Trace_QuTree_Recur(qA%Quad22, Depth+1)
      !$OMP END TASK

      !$OMP TASKWAIT
      Trace = Trace00+Trace11
    ENDIF

  END FUNCTION SpAMM_Trace_QuTree_Recur

  !> Recursive part of trace for quadtree product, @f$ \mathrm{Tr} \left[ A
  !! \times B \right] @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qB Pointer to matrix B.
  !! @param threshold The multiply threshold.
  !! @param Depth The current tier.
  !!
  !! @return The trace of the matrix produce, @f$ \mathrm{Tr} \left[ A \times B
  !! \right] @f$.
  RECURSIVE FUNCTION SpAMM_Trace_QuTree_Product_Recur(qA,qB, threshold, Depth) RESULT(Trace)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA,qB
    REAL(SpAMM_KIND), INTENT(IN)      :: threshold
    INTEGER                           :: Depth

    INTEGER :: I
    REAL(SpAMM_KIND) :: Trace
    REAL(SpAMM_KIND) :: Trace_00_00, Trace_01_10, Trace_10_01, Trace_11_11

    Trace = SpAMM_Zero

    IF(.NOT.ASSOCIATED(qA)) RETURN
    IF(.NOT.ASSOCIATED(qB)) RETURN

    IF(qA%Norm*qB%Norm < threshold) RETURN

    IF(Depth == SpAMM_TOTAL_DEPTH)THEN
      DO I = 1, SPAMM_BLOCK_SIZE
        Trace = Trace+DOT_PRODUCT(qA%Blok(I, 1:SPAMM_BLOCK_SIZE), qB%Blok(1:SPAMM_BLOCK_SIZE, I))
      ENDDO
    ELSE
      !$OMP TASK UNTIED SHARED(qA,qB,Trace_00_00)
      Trace_00_00 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad11, qB%Quad11, threshold, Depth+1)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,qB,Trace_10_01)
      Trace_10_01 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad21, qB%Quad12, threshold, Depth+1)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,qB,Trace_01_10)
      Trace_01_10 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad12, qB%Quad21, threshold, Depth+1)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,qB,Trace_11_11)
      Trace_11_11 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad22, qB%Quad22, threshold, Depth+1)
      !$OMP END TASK

      !$OMP TASKWAIT

      Trace = Trace_00_00+Trace_10_01+Trace_01_10+Trace_11_11
    ENDIF

  END FUNCTION SpAMM_Trace_QuTree_Product_Recur

  !> QuTree filter: \f$ \tilde{A} = \mathrm{filter} [A, \tau] \f$.
  RECURSIVE SUBROUTINE SpAMM_Filter_QuTree_Recur(qA,Tau,Depth)
    TYPE(QuTree), POINTER  :: qA
    REAL(SpAMM_KIND)       :: Tau
    INTEGER                :: Depth

    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(qA%Norm<Tau)THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA)
      !$OMP END TASK
      !$OMP TASKWAIT
      CALL Delete(qA)
    ELSE
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad11,Tau,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad12,Tau,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad21,Tau,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad22,Tau,Depth+1)
      !$OMP END TASK
    ENDIF
  END SUBROUTINE SpAMM_Filter_QuTree_Recur

  !> Calculate the Frobenius norm recursively on the tree.
  !!
  !! @param qA Pointer to quadtree A.
  !!
  !! @return The norm.
  RECURSIVE FUNCTION SpAMM_Norm_Reduce_QuTree_Recur(qA) RESULT(Norm)
    TYPE(QuTree), POINTER :: qA
    REAL(SpAMM_KIND)      :: Norm,Norm00,Norm01,Norm10,Norm11
    INTEGER               :: i, j

    IF(.NOT.ASSOCIATED(qA))THEN
      Norm = SpAMM_Zero
    ELSEIF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      Norm = SUM(qA%Blok**2)
      qA%Norm = SQRT(Norm)
    ELSE
      !$OMP TASK UNTIED SHARED(qA,Norm00)
      Norm00 = SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad11)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Norm01)
      Norm01 = SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad12)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Norm10)
      Norm10 = SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad21)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Norm11)
      Norm11 = SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad22)
      !$OMP END TASK
      !$OMP TASKWAIT
      Norm = Norm00+Norm01+Norm10+Norm11
      qA%Norm = SQRT(Norm)
    ENDIF
  END FUNCTION SpAMM_Norm_Reduce_QuTree_Recur

  !> Recursive linear algebra routines on row tree vectors
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree_Recur(bC,qA,bB,Depth)

    TYPE(QuTree), POINTER :: qA
    TYPE(BiTree), POINTER :: bB,bC
    INTEGER               :: Depth

    ! Associated
    IF(ASSOCIATED(qA).AND.ASSOCIATED(bB))THEN
      ! Estimate
      IF(qA%Norm*bB%Norm<SpAMM_Threshold_Multiply_QuTree_x_BiTree)RETURN
      IF(.NOT.ASSOCIATED(bC))THEN
        !$OMP CRITICAL
        ALLOCATE(bC)
        !$OMP END CRITICAL
      ENDIF
      ! Blocks
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        ! Allocate
        IF(.NOT.ALLOCATED(bC%Vect))THEN
          !$OMP CRITICAL
          ALLOCATE(bC%Vect(1:SPAMM_BLOCK_SIZE))
          bC%Vect=SpAMM_Zero
          !$OMP END CRITICAL
        END IF
        ! Accumulate
        bC%Vect(1:SPAMM_BLOCK_SIZE)=bC%Vect(1:SPAMM_BLOCK_SIZE)+MATMUL( &
          qA%Blok(1:SPAMM_BLOCK_SIZE,1:SPAMM_BLOCK_SIZE),bB%Vect(1:SPAMM_BLOCK_SIZE))
      ELSE

        ! 0=00*0
        !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(qA%Quad11%Norm*bB%Sect0%Norm>SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect0,qA%Quad11,bB%Sect0,Depth+1)
        !$OMP END TASK

        ! 1=10*0
        !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(qA%Quad21%Norm*bB%Sect0%Norm>SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect1,qA%Quad21,bB%Sect0,Depth+1)
        !$OMP END TASK

        !$OMP TASKWAIT

        ! 0=00*0+01*1
        !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(qA%Quad12%Norm*bB%Sect1%Norm>SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect0,qA%Quad12,bB%Sect1,Depth+1)
        !$OMP END TASK

        ! 1=10*0+11*1
          !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(qA%Quad22%Norm*bB%Sect1%Norm>SpAMM_RECURSION_NORMD_CUTOFF)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect1,qA%Quad22,bB%Sect1,Depth+1)
        !$OMP END TASK

        !$OMP TASKWAIT
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree_Recur

  !> BiTree in place add: \f$ A \leftarrow \alpha A + \beta B \f$.
  RECURSIVE SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA,bB,Depth)
    TYPE(BiTree),POINTER :: bA,bB
    INTEGER              :: Depth
    LOGICAL              :: TA, TB
    TA=ASSOCIATED(bA)
    TB=ASSOCIATED(bB)
    IF(TA.AND.TB)THEN
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        bA%Vect=SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha*bA%Vect &
               +SpAMM_Add_BiTree_2_BiTree_InPlace_Beta *bB%Vect
      ELSE
        !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA%Sect0,bB%Sect0,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA%Sect1,bB%Sect1,Depth+1)
        !$OMP END TASK
      ENDIF
      !$OMP TASKWAIT !! << WTF IS THIS TASKWAIT IMPORTANT FOR IFORT??
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bB,bA,Depth)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace_Recur

  !> BiTree add with replacement: \f$ C \leftarrow \alpha A + \beta B \f$.
  RECURSIVE SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA,bB,bC,Depth)

    TYPE(BiTree),POINTER :: bA,bB,bC
    INTEGER              :: Depth
    LOGICAL              :: TA, TB
    TA=ASSOCIATED(bA)
    TB=ASSOCIATED(bB)
    IF(TA.AND.TB)THEN
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        bC%Vect=SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha*bA%Vect &
          +SpAMM_Add_BiTree_2_BiTree_RePlace_Beta *bB%Vect
      ELSE
        !$OMP TASK UNTIED SHARED(bA,bB,bC) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA%Sect0,bB%Sect0,bC%Sect0,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(bA,bB,bC) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA%Sect1,bB%Sect1,bC%Sect1,Depth+1)
        !$OMP END TASK
      ENDIF
      !$OMP TASKWAIT !! << WTF IS THIS TASKWAIT IMPORTANT FOR IFORT??
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(bC,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bB,bC,Depth)
      !$OMP END TASK
    ELSEIF(.NOT.TB.AND.TA)THEN
      !$OMP TASK UNTIED SHARED(bC,bA) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA,bC,Depth)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace_Recur

  !> \f$ L_2 \f$ norm for BiTrees.
  RECURSIVE FUNCTION SpAMM_Norm_Reduce_BiTree_Recur(bA,Depth) RESULT(Norm)
    TYPE(BiTree), POINTER :: bA
    INTEGER               :: Depth
    REAL(SpAMM_KIND)      :: Norm, Norm0, Norm1

    IF(.NOT.ASSOCIATED(bA))THEN
       Norm=SpAMM_Zero
       RETURN
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
       Norm=SUM(bA%Vect(1:SPAMM_BLOCK_SIZE)**2)
       bA%Norm=SQRT(Norm)
    ELSE
       !$OMP TASK UNTIED SHARED(bA) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Norm0=SpAMM_Norm_Reduce_BiTree_Recur(bA%Sect0,Depth+1)
       !$OMP END TASK
       !$OMP TASK UNTIED SHARED(bA) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Norm1=SpAMM_Norm_Reduce_BiTree_Recur(bA%Sect1,Depth+1)
       !$OMP END TASK
       !$OMP TASKWAIT
       Norm=Norm0+Norm1
       bA%Norm=SQRT(Norm)
    ENDIF
  END FUNCTION SpAMM_Norm_Reduce_BiTree_Recur

  !> Dot product for BiTrees
  RECURSIVE FUNCTION SpAMM_Dot_Product_BiTree_Recur(bA,bB,Depth) RESULT(Dot)
    TYPE(BiTree), POINTER :: bA,bB
    INTEGER               :: Depth
    REAL(SpAMM_KIND)      :: Dot, Dot0, Dot1

    IF(.NOT.ASSOCIATED(bA))THEN
       Dot=SpAMM_Zero
       RETURN
    ELSEIF(.NOT.ASSOCIATED(bB))THEN
       Dot=SpAMM_Zero
       RETURN
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
       Dot=DOT_PRODUCT(bA%Vect(1:SPAMM_BLOCK_SIZE),bB%Vect(1:SPAMM_BLOCK_SIZE))
    ELSE
       !$OMP TASK UNTIED SHARED(bA,bB,Dot0) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Dot0=SpAMM_Dot_Product_BiTree_Recur(bA%Sect0,bB%Sect0,Depth+1)
       !$OMP END TASK
       !$OMP TASK UNTIED SHARED(bA,bB,Dot1) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Dot1=SpAMM_Dot_Product_BiTree_Recur(bA%Sect1,bB%Sect1,Depth+1)
       !$OMP END TASK
       !$OMP TASKWAIT
       Dot=Dot0+Dot1
    ENDIF
  END FUNCTION SpAMM_Dot_Product_BiTree_Recur

  !> Multiply vector with scalar.
  RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar_Recur(bA,a,Depth)
    INTEGER :: Depth
    TYPE(BiTree), POINTER   :: bA
    REAL(SpAMM_KIND)        :: a
    IF(.NOT.ASSOCIATED(bA))RETURN
    IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(bA%Vect))THEN
       bA%Norm=bA%Norm*ABS(a)
       bA%Vect=bA%Vect*a
    ELSE
       !$OMP TASK UNTIED SHARED(bA) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA%Sect0,a,Depth+1)
       !$OMP END TASK
       !$OMP TASK UNTIED  SHARED(bA) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA%Sect1,a,Depth+1)
       !$OMP END TASK
       bA%Norm=bA%Norm*ABS(a)
    ENDIF
  END SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar_Recur

  !> Multiply two 2nd order matrices: \f$ C \leftarrow \alpha A B + \beta C \f$.
  !!
  !! If the tolerance is not given, then \f$ \tau = 0 \f$ is used, i.e. the product reverts to an exact dense product.
  !!
  !! @bug The code for \f$ \alpha \neq 1 \f$ and \f$ \beta \neq 0 \f$ is not implemented yet.
  !!
  !! @param A The matrix \f$ A \f$.
  !! @param B The matrix \f$ B \f$.
  !! @param C The matrix \f$ C \f$.
  !! @param tolerance The SpAMM tolerance \f$ \tau \f$.
  !! @param alpha The scalar \f$ \alpha \f$.
  !! @param beta The scalar \f$ \beta \f$.
  subroutine spamm_multiply_2nd_order_x_2nd_order (A, B, C, tolerance, alpha, beta)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A, B
    type(spamm_matrix_2nd_order), pointer, intent(inout) :: C
    real(spamm_kind), intent(in), optional :: tolerance
    real(spamm_kind), intent(in), optional :: alpha, beta

    real(spamm_kind) :: local_tolerance

    if(present(tolerance)) then
      local_tolerance = tolerance
    else
      local_tolerance = 0
    endif

    call spamm_multiply_qutree_x_qutree(A%root, B%root, C%root, tolerance)

  end subroutine spamm_multiply_2nd_order_x_2nd_order

  !> Frobenius norm of 2nd order matrix. This function updates the norm on the matrix and returns the square of the norm.
  !!
  !! @param A The matrix.
  !!
  !! @return The squared Frobenius norm.
  function spamm_norm_reduce_matrix_2nd_order (A) result (norm)

    real(spamm_kind) :: norm
    type(spamm_matrix_2nd_order), pointer, intent(in) :: A

    norm = spamm_norm_reduce_qutree_recur(A%root)

  end function spamm_norm_reduce_matrix_2nd_order

END MODULE SpAMM_ALGEBRA
