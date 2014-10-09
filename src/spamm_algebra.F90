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
!! @author Nicolas Bock nicolasbock@freeon.org
MODULE SpAMM_ALGEBRA

  use spamm_types
  use spamm_globals
  use spamm_management
  use spamm_utilities

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
  REAL(SpAMM_KIND) :: SpAMM_Threshold_Multiply_QuTree_x_BiTree

  !> Interface for multiplication operations between different SpAMM types.
  !!
  !! The Sparse Approximate Matrix-Multiply (SpAMM):
  !! @f$ C \leftarrow A \cdot B @f$.
  INTERFACE Multiply
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_QuTree
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_Scalar
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_BiTree
    MODULE PROCEDURE SpAMM_Multiply_BiTree_x_Scalar
    module procedure spamm_multiply_order_1_x_scalar
    module procedure spamm_multiply_2nd_order_x_2nd_order
    module procedure spamm_multiply_2nd_order_x_scalar
    module procedure spamm_multiply_2nd_order_x_1st_order
  END INTERFACE

  !> Interface for trace operations.
  INTERFACE Trace
    MODULE PROCEDURE SpAMM_Trace_QuTree
    MODULE PROCEDURE SpAMM_Trace_QuTree_Product
    module procedure spamm_trace_2nd_order
    module procedure spamm_trace_2nd_order_product
  END INTERFACE

  !> Interface for additions operations between different SpAMM types.
  INTERFACE Add
    MODULE PROCEDURE SpAMM_Add_QuTree_2_QuTree_InPlace
    MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_InPlace
    MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_RePlace
    MODULE PROCEDURE SpAMM_Add_Identity_2_QuTree_InPlace
    module procedure spamm_add_identity_to_matrix_2nd_order
    module procedure spamm_add_2nd_order_to_2nd_order
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

  !> Multiplication operation between a quadtree and a quadtree, @f$ C \leftarrow A \cdot B @f$.
  !!
  !! @param qA Pointer to quadtree A.
  !! @param qB Pointer to quadtree B.
  !! @param qC Pointer to quadtree C.
  !! @param threshold The SpAMM threshold.
  SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree(qA, qB, qC, threshold)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA, qB
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
    REAL(SpAMM_KIND), OPTIONAL :: threshold

    real(spamm_kind) :: local_threshold
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()

    if(.not. associated(qA) .or. .not. associated(qB)) then
      LOG_DEBUG("either A or B are not associated")
      return
    endif

    if(.not. associated(qC)) then
      CALL NewQuNode(qC, qA%i_lower, qB%j_lower, qA%i_upper, qB%j_upper)
    endif

    !$OMP PARALLEL

    ! The master thread will lead execution of the product. All subsequent tasks
    ! are untied and can be executed by any thread in the thread group.
    !$OMP MASTER

#ifdef _OPENMP
    LOG_INFO("Multiply on "//to_string(omp_get_num_threads())//" OpenMP threads")
#endif

    IF(PRESENT(threshold))THEN
      local_threshold = threshold
    ELSE
      local_threshold = 0
    ENDIF

    qC%number_operations = 0

    LOG_DEBUG("resetting C")
    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, SpAMM_Zero)

    LOG_DEBUG("recursive multiplication")

    !$OMP TASK UNTIED SHARED(qA,qB,qC)
    CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC, qA, qB, local_threshold)
    !$OMP END TASK

    !$OMP END MASTER

    !$OMP END PARALLEL

    qC%norm = norm(qC)
    qC%norm = sqrt(qC%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree

  !> Scalar multiply of 2nd order matrix: @f$ A \leftarrow \alpha A @f$.
  !!
  !! @param A The matrix.
  !! @param alpha The scalar \f$ \alpha \f$.
  subroutine spamm_multiply_2nd_order_x_scalar (A, alpha)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A
    real(spamm_kind), intent(in) :: alpha

    LOG_DEBUG("multiplying matrix by scalar "//to_string(alpha))
    call spamm_multiply_qutree_x_scalar(A%root, alpha)
    A%number_operations = A%root%number_operations

  end subroutine spamm_multiply_2nd_order_x_scalar

  !> Scalar multiply: @f$ A \leftarrow alpha A @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param alpha Scalar @f$ \alpha @f$.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar(qA, alpha)

    TYPE(QuTree), POINTER        :: qA
    REAL(SpAMM_KIND), INTENT(IN) :: alpha

    INTEGER            :: Depth
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    IF(.NOT.ASSOCIATED(qA)) then
      LOG_DEBUG("qA not associated")
      RETURN
    endif

    Depth=0
    TInitial = SpAMM_Get_Time()

    qA%number_operations = 0

    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)
    !$OMP END TASK

    !$OMP TASKWAIT

    qA%norm = norm(qA)
    qA%norm = sqrt(qA%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_Scalar",3)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar

  !> Add 2 matrices, @f$ A \leftarrow \alpha A + \beta B @f$.
  !!
  !! If any of the two factors Alpha or Beta are not supplied then they default
  !! to one.
  !!
  !! @param A Pointer to matrix A.
  !! @param B Pointer to matrix B.
  !! @param alpha Factor @f$ \alpha @f$.
  !! @param beta Factor @f$ \beta @f$.
  subroutine spamm_add_2nd_order_to_2nd_order (A, B, alpha, beta)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A
    type(spamm_matrix_2nd_order), pointer, intent(in) :: B
    real(spamm_kind), intent(in), optional :: alpha, beta

    LOG_DEBUG("Adding matrices: alpha = "//to_string(alpha)//", beta = "//to_string(beta))

    if(.not. associated(B)) then
      return
    endif

    if(.not. associated(A)) then
      call spamm_allocate_matrix_2nd_order(B%M, B%N, A)
    endif

    call spamm_add_qutree_2_qutree_inplace_recur(A%root, B%root, alpha, beta)

    A%norm = A%root%norm
    A%number_nonzeros = A%root%number_nonzeros

  end subroutine spamm_add_2nd_order_to_2nd_order

  !> Add 2 quadtree matrices, @f$ A \leftarrow \alpha A + \beta B @f$.
  !!
  !! If any of the two factors Alpha or Beta are not supplied then they default
  !! to one.
  !!
  !! @param qA [inout] Pointer to matrix A.
  !! @param qB [in] Pointer to matrix B.
  !! @param alpha Factor @f$ \alpha @f$.
  !! @param beta Factor @f$ \beta @f$.
  SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace (qA, qB, Alpha, Beta)

    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    TYPE(QuTree), POINTER, INTENT(IN) :: qB
    REAL(SpAMM_KIND), intent(in), OPTIONAL :: Alpha, Beta
    REAL(SpAMM_KIND) :: InPlace_Alpha, InPlace_Beta
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

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

    LOG_DEBUG("Adding "//to_string(inplace_alpha)//"*A + "//to_string(inplace_beta)//"*B")

    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(qA,qB)
    CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA, qB, InPlace_Alpha, InPlace_Beta)
    !$OMP END TASK
    !$OMP TASKWAIT

    qA%norm = norm(qA)
    qA%norm = sqrt(qA%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_QuTree_2_QuTree_InPlace",4)

  END SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace

  !> In-place add for 2nd order SpAMM matrix: @f$ A \leftarrow \alpha \mathrm{Id} @f$.
  !!
  !! @param alpha The factor \f$ \alpha \f$.
  subroutine spamm_add_identity_to_matrix_2nd_order (A, alpha)

    type(spamm_matrix_2nd_order), pointer, intent(inout) :: A
    real(spamm_kind), intent(in) :: alpha

    call spamm_add_identity_2_qutree_inplace_recur(A%root, alpha, A%M, A%N)

  end subroutine spamm_add_identity_to_matrix_2nd_order

  !> QuTree In Place Add: \f$ A \leftarrow A + \alpha \mathrm{Id} \f$.
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

  !> Trace for QuTree: @f$ a = \mathrm{Tr} [ A ] @f$.
  !!
  !! @param qA Pointer to matrix A.
  !!
  !! @return The trace.
  FUNCTION SpAMM_Trace_QuTree (qA) RESULT(a)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA
    REAL(SpAMM_KIND)                  :: a

    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()
    !$OMP TASK UNTIED SHARED(qA,a)
    a = SpAMM_Trace_QuTree_Recur(qA)
    !$OMP END TASK
    !$OMP TASKWAIT

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree",6)

  END FUNCTION SpAMM_Trace_QuTree

  !> The trace of a 2nd order matrix: @f$ a \leftarrow \mathrm{Tr} [ A ] @f$.
  !!
  !! @param A The matrix.
  !!
  !! @return The trace.
  function spamm_trace_2nd_order (A) result (trace_a)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    real(spamm_kind) :: trace_a

    trace_a = 0

    if(.not. associated(A)) then
      return
    endif

    if(.not. associated(A%root)) then
      return
    endif

    trace_a = spamm_trace_qutree(A%root)

  end function spamm_trace_2nd_order

  !> Trace for quadtree product, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qB Pointer to matrix B.
  !!
  !! @return The trace of the matrix produce, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
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
      multiplyThreshold = 0
    ENDIF

    TInitial = SpAMM_Get_Time()

    !$OMP TASK UNTIED SHARED(qA,a)
    a = SpAMM_Trace_QuTree_Product_Recur(qA, qB, multiplyThreshold)
    !$OMP END TASK
    !$OMP TASKWAIT

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree_Product",7)

  END FUNCTION SpAMM_Trace_QuTree_Product

  !> Trace for matrix product, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
  !!
  !! @param A Pointer to matrix A.
  !! @param B Pointer to matrix B.
  !!
  !! @return The trace of the matrix produce, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
  function spamm_trace_2nd_order_product (A, B, tolerance) result (trace_ab)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A, B
    real(spamm_kind), intent(in), optional :: tolerance
    real(spamm_kind) :: trace_ab

    real(spamm_kind) :: local_tolerance

    trace_ab = 0

    if(.not. associated(A) .or. .not. associated(B)) then
      LOG_INFO("either A or B are not associated")
      return
    endif

    if(.not. associated(A%root) .or. .not. associated(B%root)) then
      LOG_INFO("either A%root or B%root are not associated")
      return
    endif

    if(present(tolerance)) then
      local_tolerance = tolerance
    else
      local_tolerance = 0
    endif

    trace_ab = spamm_trace_qutree_product(A%root, B%root, local_tolerance)

  end function spamm_trace_2nd_order_product

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
    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Norm_Reduce_QuTree",9)

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

  !> Sparse Approximate Matrix-Vector Multiply (SpAMV): \f$ C \leftarrow A \cdot B \f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param bB Pointer to vector B.
  !! @param bC Pointer to vector C.
  !! @param tolerance The SpAMM tolerance.
  SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree (qA, bB, bC, tolerance)

    TYPE(QuTree), POINTER, intent(in) :: qA
    TYPE(BiTree), POINTER, intent(in) :: bB
    TYPE(BiTree), POINTER, intent(inout) :: bC
    REAL(SpAMM_KIND),OPTIONAL :: tolerance
    real(spamm_kind) :: local_tolerance
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    IF(PRESENT(tolerance))THEN
      local_tolerance = tolerance
    ELSE
      local_tolerance = 0
    ENDIF

    TInitial = SpAMM_Get_Time()

    CALL SpAMM_Multiply_BiTree_x_Scalar(bC, SpAMM_Zero)
    !$OMP TASK UNTIED SHARED(qA,bB,bC)
    CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC, qA, bB, local_tolerance)
    !$OMP END TASK
    !$OMP TASKWAIT

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_BiTree",13)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree

  !> Scalar multiply: \f$ A \leftarrow \alpha A \f$.
  !!
  !! @param bA Pointer to vector A.
  !! @param alpha Scalar @f$ \alpha @f$.
  RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar(bA, alpha)

    TYPE(BiTree), POINTER :: bA
    REAL(SpAMM_KIND) :: alpha
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    IF(.NOT.ASSOCIATED(bA))RETURN

    TInitial = SpAMM_Get_Time()
    !$OMP TASK SHARED(bA)
    CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA,alpha)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_BiTree_x_Scalar",14)

  END SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar

  !> Multiply a vector with a scalar: @f$ V \leftarrow \alpha V @f$.
  !!
  !! @param V The vector
  !! @param alpha The scalar @f$ \alpha @f$
  subroutine spamm_multiply_order_1_x_scalar (V, alpha)

    type(spamm_matrix_order_1), pointer, intent(inout) :: V
    real(spamm_kind), intent(in) :: alpha

    call spamm_multiply_bitree_x_scalar(V%root, alpha)

  end subroutine spamm_multiply_order_1_x_scalar

  !> Norm for BiTrees.
  !!
  !! @param bA Pointer to vector A.
  !!
  !! @result The norm of the vector.
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
  !! quadtree, @f$ C \leftarrow A \cdot B @f$.
  !!
  !! @param qA Pointer to quadtree A.
  !! @param qB Pointer to quadtree B.
  !! @param qC Pointer to quadtree C.
  !! @param threshold The SpAMM product tolerance.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur(qC, qA, qB, threshold)

    TYPE(QuTree), POINTER :: qC, qA, qB
    REAL(SpAMM_KIND) :: threshold
    !real(spamm_kind) :: temp
    integer :: i, j

    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB)) THEN
      LOG_DEBUG("qA: "//to_string(qA))
      LOG_DEBUG("qB: "//to_string(qB))

      ! Apply the SpAMM condition.
      if(qA%Norm*qB%Norm <= threshold) then
        LOG_DEBUG("going back up")
        RETURN
      endif

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

      LOG_DEBUG("qC: "//to_string(qC))
      LOG_DEBUG("   operations = "//to_string(qC%number_operations))

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
#ifdef SPAMM_COUNTERS
        qC%number_operations = qC%number_operations+SPAMM_BLOCK_SIZE**3
        qC%number_nonzeros = sum(reshape( &
          [ ((1, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) ], &
          [ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE ]), &
          reshape( &
          [ ((qC%blok(i, j) /= 0.0, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) ], &
          [ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE ]))
#endif

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

#ifdef SPAMM_COUNTERS
        qC%number_operations = 0
        qC%number_nonzeros = 0

        if(associated(qC%quad11)) then
          qC%number_operations = qC%number_operations+qC%Quad11%number_operations
          qC%number_nonzeros = qC%number_nonzeros+qC%Quad11%number_nonzeros
        endif

        if(associated(qC%quad12)) then
          qC%number_operations = qC%number_operations+qC%Quad12%number_operations
          qC%number_nonzeros = qC%number_nonzeros+qC%Quad12%number_nonzeros
        endif

        if(associated(qC%quad21)) then
          qC%number_operations = qC%number_operations+qC%Quad21%number_operations
          qC%number_nonzeros = qC%number_nonzeros+qC%Quad21%number_nonzeros
        endif

        if(associated(qC%quad22)) then
          qC%number_operations = qC%number_operations+qC%Quad22%number_operations
          qC%number_nonzeros = qC%number_nonzeros+qC%Quad22%number_nonzeros
        endif
#endif

      ENDIF
    else
      LOG_DEBUG("either A or B is not associated")
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur

  !> Recursive part of scalar multiply with quadtree matrix, @f$ A \leftarrow \alpha A @f$.
  !!
  !! @param qA Pointer to quadtree.
  !! @param alpha The scalar.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)

    TYPE(QuTree), POINTER :: qA
    REAL(SpAMM_KIND) :: alpha

    IF(.NOT.ASSOCIATED(qA)) then
      LOG_DEBUG("qA not associated")
      RETURN
    endif

    LOG_DEBUG("qA = "//to_string(qA))

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      ! At the bottom, multiply the block.
      qA%Norm = qA%Norm*ABS(alpha)
      qA%Blok = alpha*qA%Blok
#ifdef SPAMM_COUNTERS
      qA%number_operations = SPAMM_BLOCK_SIZE**2
#endif
    ELSE
      !$OMP TASK UNTIED SHARED(qA)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad11, alpha)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA)
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

#ifdef SPAMM_COUNTERS
      qA%number_operations = 0
      qA%number_nonzeros = 0

      if(associated(qA%quad11)) then
        qA%number_operations = qA%number_operations+qA%Quad11%number_operations
        qA%number_nonzeros = qA%number_nonzeros+qA%Quad11%number_nonzeros
      endif

      if(associated(qA%quad12)) then
        qA%number_operations = qA%number_operations+qA%Quad12%number_operations
        qA%number_nonzeros = qA%number_nonzeros+qA%Quad12%number_nonzeros
      endif

      if(associated(qA%quad21)) then
        qA%number_operations = qA%number_operations+qA%Quad21%number_operations
        qA%number_nonzeros = qA%number_nonzeros+qA%Quad21%number_nonzeros
      endif

      if(associated(qA%quad22)) then
        qA%number_operations = qA%number_operations+qA%Quad22%number_operations
        qA%number_nonzeros = qA%number_nonzeros+qA%Quad22%number_nonzeros
      endif
#endif
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur

  !> Add 2 quadtree matrices, @f$ A \leftarrow \alpha A + \beta B @f$.
  !!
  !! @param qA [inout] Pointer to matrix A.
  !! @param qB [in] Pointer to matrix B.
  !! @param alpha Factor @f$ \alpha @f$.
  !! @param beta Factor @f$ \beta @f$.
  RECURSIVE SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA, qB, alpha, beta)

    TYPE(QuTree), POINTER :: qA,qB
    REAL(SpAMM_KIND) :: alpha, beta
    LOGICAL :: TA, TB
    integer :: i, j

    TA=ASSOCIATED(qA)
    TB=ASSOCIATED(qB)

    IF(TA.AND.TB)THEN
      LOG_DEBUG("qA = "//to_string(qA))
      LOG_DEBUG("qB = "//to_string(qB))

      IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) THEN
        qA%Blok = alpha*qA%Blok+beta*qB%Blok
        qA%norm = sqrt(sum(qA%blok**2))
#ifdef SPAMM_COUNTERS
        qA%number_nonzeros = sum(reshape( &
          [ ((1, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) ], &
          [ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE ]), &
          reshape( &
          [ ((qA%blok(i, j) /= 0.0, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) ], &
          [ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE ]))
#endif
      ELSE
        !$OMP TASK UNTIED SHARED(qA,qB)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad11, qB%Quad11, alpha, beta)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad12, qB%Quad12, alpha, beta)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad21, qB%Quad21, alpha, beta)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad22, qB%Quad22, alpha, beta)
        !$OMP END TASK
        !$OMP TASKWAIT

        qA%norm = 0
        qA%number_nonzeros = 0

        if(associated(qA%quad11)) then
          qA%norm = qA%norm+qA%quad11%norm**2
          qA%number_nonzeros = qA%number_nonzeros+qA%quad11%number_nonzeros
        endif

        if(associated(qA%quad12)) then
          qA%norm = qA%norm+qA%quad12%norm**2
          qA%number_nonzeros = qA%number_nonzeros+qA%quad12%number_nonzeros
        endif

        if(associated(qA%quad21)) then
          qA%norm = qA%norm+qA%quad21%norm**2
          qA%number_nonzeros = qA%number_nonzeros+qA%quad21%number_nonzeros
        endif

        if(associated(qA%quad22)) then
          qA%norm = qA%norm+qA%quad22%norm**2
          qA%number_nonzeros = qA%number_nonzeros+qA%quad22%number_nonzeros
        endif

        qA%norm = sqrt(qA%norm)
      ENDIF
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(qA,qB)
      CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qB, qA)
      call spamm_multiply_qutree_x_scalar_recur(qA, beta)
      !$OMP END TASK
      !$OMP TASKWAIT
    ELSEIF(TA .AND. .NOT.TB) THEN
      ! Multiply A tree with alpha.
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)
    else
      LOG_DEBUG("either A or B are not associated")
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
    integer :: i
    !integer :: j

    LOG_DEBUG("q:"//to_string(qA%i_lower)//" "//to_string(qA%j_lower))
    LOG_DEBUG("  "//to_string(qA%i_upper)//" "//to_string(qA%j_upper))

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
      if(.not. allocated(qA%blok)) then
        ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
        qA%Blok = SpAMM_Zero
      ENDIF

      !write(*, *) "before"
      !do i = 1, SPAMM_BLOCK_SIZE
      !  write(*, "(4f10.3)") (qA%blok(i, j), j = 1, SPAMM_BLOCK_SIZE)
      !enddo

      do i = 1, MIN(SPAMM_BLOCK_SIZE, M-qA%i_lower+1, N-qA%j_lower+1)
        qA%Blok(i:i, i:i) = qA%Blok(i:i, i:i)+alpha
      enddo

      !write(*, *) "after"
      !do i = 1, SPAMM_BLOCK_SIZE
      !  write(*, "(4f10.3)") (qA%blok(i, j), j = 1, SPAMM_BLOCK_SIZE)
      !enddo
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
  !!
  !! @return The trace.
  RECURSIVE FUNCTION SpAMM_Trace_QuTree_Recur(qA) RESULT(Trace)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA
    REAL(SpAMM_KIND)                  :: Trace

    REAL(SpAMM_KIND) :: Trace00, Trace11
    INTEGER          :: I

    Trace = SpAMM_Zero

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) THEN
      Trace = SpAMM_Zero
      IF(.NOT.ASSOCIATED(qA)) RETURN
      DO I = 1, SPAMM_BLOCK_SIZE
        Trace = Trace+qA%Blok(I,I)
      ENDDO
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11).AND. &
        .NOT.ASSOCIATED(qA%Quad22))THEN
      Trace = SpAMM_Zero
    ELSEIF(.NOT.ASSOCIATED(qA%Quad22))THEN
      !$OMP TASK UNTIED SHARED(qA)
      Trace = SpAMM_Trace_QuTree_Recur(qA%Quad11)
      !$OMP END TASK
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11))THEN
      !$OMP TASK UNTIED SHARED(qA)
      Trace = SpAMM_Trace_QuTree_Recur(qA%Quad22)
      !$OMP END TASK
    ELSE
      !$OMP TASK UNTIED SHARED(qA,Trace00)
      Trace00 = SpAMM_Trace_QuTree_Recur(qA%Quad11)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,Trace11)
      Trace11 = SpAMM_Trace_QuTree_Recur(qA%Quad22)
      !$OMP END TASK

      !$OMP TASKWAIT
      Trace = Trace00+Trace11
    ENDIF

  END FUNCTION SpAMM_Trace_QuTree_Recur

  !> Recursive part of trace for quadtree product, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
  !!
  !! @param qA Pointer to matrix A.
  !! @param qB Pointer to matrix B.
  !! @param threshold The multiply threshold.
  !!
  !! @return The trace of the matrix produce, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
  RECURSIVE FUNCTION SpAMM_Trace_QuTree_Product_Recur(qA,qB, threshold) RESULT(Trace)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA,qB
    REAL(SpAMM_KIND), INTENT(IN)      :: threshold

    INTEGER :: I
    REAL(SpAMM_KIND) :: Trace
    REAL(SpAMM_KIND) :: Trace_00_00, Trace_01_10, Trace_10_01, Trace_11_11

    Trace = SpAMM_Zero

    IF(.NOT.ASSOCIATED(qA)) RETURN
    IF(.NOT.ASSOCIATED(qB)) RETURN

    IF(qA%Norm*qB%Norm < threshold) RETURN

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE)THEN
      DO I = 1, SPAMM_BLOCK_SIZE
        Trace = Trace+DOT_PRODUCT(qA%Blok(I, 1:SPAMM_BLOCK_SIZE), qB%Blok(1:SPAMM_BLOCK_SIZE, I))
      ENDDO
    ELSE
      !$OMP TASK UNTIED SHARED(qA,qB,Trace_00_00)
      Trace_00_00 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad11, qB%Quad11, threshold)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,qB,Trace_10_01)
      Trace_10_01 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad21, qB%Quad12, threshold)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,qB,Trace_01_10)
      Trace_01_10 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad12, qB%Quad21, threshold)
      !$OMP END TASK

      !$OMP TASK UNTIED SHARED(qA,qB,Trace_11_11)
      Trace_11_11 = SpAMM_Trace_QuTree_Product_Recur(qA%Quad22, qB%Quad22, threshold)
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

  !> Recursive linear algebra routines on row tree vectors: @f$ C \leftarrow A \dot B @f$.
  !!
  !! @param bC Vector C.
  !! @param qA Matrix A.
  !! @param bB Vector B.
  !! @param tolerance The SpAMM tolerance.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree_Recur(bC, qA, bB, tolerance)

    TYPE(QuTree), POINTER, intent(in) :: qA
    TYPE(BiTree), POINTER, intent(in) :: bB
    TYPE(BiTree), POINTER, intent(inout) :: bC
    real(spamm_kind), intent(in) :: tolerance

    ! Associated
    IF(ASSOCIATED(qA).AND.ASSOCIATED(bB)) THEN
      ! Estimate
      LOG_DEBUG("norm(A) = "//to_string(qA%norm))
      LOG_DEBUG("norm(B) = "//to_string(bB%norm))

      IF(qA%Norm*bB%Norm < tolerance) RETURN

      IF(.NOT.ASSOCIATED(bC))THEN
        LOG_DEBUG("allocating new node in C bitree")
        ALLOCATE(bC)
      ENDIF

      ! Blocks
      IF(bC%i_upper-bC%i_lower+1 == SPAMM_BLOCK_SIZE) THEN
        ! Allocate
        IF(.NOT.ALLOCATED(bC%Vect))THEN
          !$OMP CRITICAL
          ALLOCATE(bC%Vect(SPAMM_BLOCK_SIZE))
          bC%Vect=SpAMM_Zero
          !$OMP END CRITICAL
        END IF
        ! Accumulate
        bC%Vect(1:SPAMM_BLOCK_SIZE)=bC%Vect+MATMUL(qA%Blok, bB%Vect)
      ELSE

        ! 0=00*0
        !$OMP TASK UNTIED SHARED(qA,bB,bC)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%sect1, qA%Quad11, bB%sect1, tolerance)
        !$OMP END TASK

        ! 1=10*0
        !$OMP TASK UNTIED SHARED(qA,bB,bC)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%sect2, qA%Quad21, bB%sect1, tolerance)
        !$OMP END TASK

        ! 0=00*0+01*1
        !$OMP TASK UNTIED SHARED(qA,bB,bC)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%sect1, qA%Quad12, bB%sect2, tolerance)
        !$OMP END TASK

        ! 1=10*0+11*1
        !$OMP TASK UNTIED SHARED(qA,bB,bC)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%sect2, qA%Quad22, bB%sect2, tolerance)
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
        CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA%sect1,bB%sect1,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA%sect2,bB%sect2,Depth+1)
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
        CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA%sect1,bB%sect1,bC%sect1,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(bA,bB,bC) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA%sect2,bB%sect2,bC%sect2,Depth+1)
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
      Norm0 = 0
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm0=SpAMM_Norm_Reduce_BiTree_Recur(bA%sect1,Depth+1)
      !$OMP END TASK

      norm1 = 0
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm1=SpAMM_Norm_Reduce_BiTree_Recur(bA%sect2,Depth+1)
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
       Dot0=SpAMM_Dot_Product_BiTree_Recur(bA%sect1,bB%sect1,Depth+1)
       !$OMP END TASK
       !$OMP TASK UNTIED SHARED(bA,bB,Dot1) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Dot1=SpAMM_Dot_Product_BiTree_Recur(bA%sect2,bB%sect2,Depth+1)
       !$OMP END TASK
       !$OMP TASKWAIT
       Dot=Dot0+Dot1
    ENDIF
  END FUNCTION SpAMM_Dot_Product_BiTree_Recur

  !> Multiply vector with scalar.
  RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar_Recur(bA, alpha)

    TYPE(BiTree), POINTER :: bA
    REAL(SpAMM_KIND) :: alpha

    IF(.NOT.ASSOCIATED(bA))RETURN

    LOG_DEBUG("q: "//to_string(bA%i_lower)//", "//to_string(bA%i_upper))

    IF(bA%i_upper-bA%i_lower+1 == SPAMM_BLOCK_SIZE .AND. ALLOCATED(bA%Vect))THEN
       bA%Norm=bA%Norm*ABS(alpha)
       bA%Vect=bA%Vect*alpha
       LOG_DEBUG("vect = "//to_string(bA%vect(1)))
    ELSE
       if(associated(bA%sect1)) then
         !$OMP TASK UNTIED SHARED(bA)
         CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA%sect1,alpha)
         !$OMP END TASK
       endif

       if(associated(bA%sect2)) then
         !$OMP TASK UNTIED  SHARED(bA)
         CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA%sect2,alpha)
         !$OMP END TASK
       endif

       bA%Norm=bA%Norm*ABS(alpha)
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
    real(spamm_kind) :: local_alpha, local_beta

    if(present(tolerance)) then
      local_tolerance = tolerance
    else
      local_tolerance = 0
    endif

    if(present(alpha)) then
      local_alpha = alpha
    else
      local_alpha = 1
    endif

    if(present(beta)) then
      local_beta = beta
    else
      local_beta = 0
    endif

    if(.not. associated(A) .or. .not. associated(B)) then
      LOG_INFO("either A or B are not allocated")
      return
    endif

    if(associated(A, C)) then
      LOG_FATAL("A == C")
      error stop
    endif

    if(associated(B, C)) then
      LOG_FATAL("B == C")
      error stop
    endif

    if(.not. associated(C)) then
      call new(A%M, B%N, C)
    endif

    LOG_DEBUG("multiplying A*B with tolerance "//to_string(local_tolerance))

    call reset_counters(C)
    call spamm_multiply_qutree_x_qutree(A%root, B%root, C%root, local_tolerance)

    if(associated(C%root)) then
      C%norm = C%root%norm
      C%number_operations = C%root%number_operations
      C%number_nonzeros = C%root%number_nonzeros
    endif

  end subroutine spamm_multiply_2nd_order_x_2nd_order

  subroutine spamm_multiply_2nd_order_x_1st_order (A, B, C, threshold)

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    type(spamm_matrix_order_1), pointer, intent(in) :: B
    type(spamm_matrix_order_1), pointer, intent(out) :: C
    real(spamm_kind), optional, intent(in) :: threshold

    real(spamm_kind) :: local_threshold

    if(.not. associated(C)) then
      call new(A%M, C)
    endif

    if(A%N /= B%N) then
      call write_log(FATAL, "dimension mismatch: A%N ("//to_string(A%N)//" /= "//"B%N ("//to_string(B%N)//")")
    endif

    if(A%M /= C%N) then
      call write_log(FATAL, "dimension mismatch: A%M ("//to_string(A%M)//" /= "//"C%N ("//to_string(C%N)//")")
    endif

    if(present(threshold)) then
      local_threshold = threshold
    else
      local_threshold = 0
    endif

    LOG_DEBUG("multiply matrix*vector")
    call spamm_multiply_qutree_x_bitree(A%root, B%root, C%root, local_threshold)

  end subroutine spamm_multiply_2nd_order_x_1st_order

  !> Frobenius norm of 2nd order matrix. This function updates the norm on the matrix and returns the square of the norm.
  !!
  !! @param A The matrix.
  !!
  !! @return The squared Frobenius norm.
  function spamm_norm_reduce_matrix_2nd_order (A) result (norm)

    real(spamm_kind) :: norm
    type(spamm_matrix_2nd_order), pointer, intent(in) :: A

    A%norm = spamm_norm_reduce_qutree_recur(A%root)
    norm = A%norm

  end function spamm_norm_reduce_matrix_2nd_order

END MODULE SpAMM_ALGEBRA
