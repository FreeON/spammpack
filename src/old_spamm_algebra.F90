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
    REAL(kind(0d0)) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()

    if(.not. associated(qA) .or. .not. associated(qB)) then
!       LOG_DEBUG("either A or B are not associated")
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
!    LOG_INFO("Multiply on "//to_string(omp_get_num_threads())//" OpenMP threads")
#endif

    IF(PRESENT(threshold))THEN
       local_threshold = threshold
    ELSE
       local_threshold = 0
    ENDIF

    qC%number_operations = 0

!    LOG_DEBUG("resetting C")
    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, 0d0)

!    LOG_DEBUG("recursive multiplication")

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

  SUBROUTINE SpAMM_Multiply_QuTreeT_x_QuTree(qA, qB, qC, threshold)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA, qB
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
    REAL(SpAMM_KIND), OPTIONAL :: threshold

    real(spamm_kind) :: local_threshold
    real(kind(0d0)) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()

    if(.not. associated(qA) .or. .not. associated(qB)) then
!       LOG_DEBUG("either A or B are not associated")
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
!    LOG_INFO("Multiply on "//to_string(omp_get_num_threads())//" OpenMP threads")
#endif

    IF(PRESENT(threshold))THEN
       local_threshold = threshold
    ELSE
       local_threshold = 0
    ENDIF

    qC%number_operations = 0

!    LOG_DEBUG("resetting C")
    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, 0d0)

!    LOG_DEBUG("recursive multiplication")

    !$OMP TASK UNTIED SHARED(qA,qB,qC)
    !    CALL SpAMM_Multiply_QuTreeT_x_QuTree_Recur(qC, qA, qB, local_threshold)
    !$OMP END TASK

    !$OMP END MASTER

    !$OMP END PARALLEL

    qC%norm = norm(qC)
    qC%norm = sqrt(qC%norm)

    TTotal=SpAMM_Get_Time()-TInitial
    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)

  END SUBROUTINE SpAMM_Multiply_QuTreeT_x_QuTree

  !> Scalar multiply of 2nd order matrix: @f$ A \leftarrow \alpha A @f$.
  !!
  !! @param A The matrix.
  !! @param alpha The scalar \f$ \alpha \f$.
  subroutine spamm_multiply_2nd_order_x_scalar (A, alpha)

    type(spamm_matrix_order_2), pointer, intent(inout) :: A
    real(spamm_kind), intent(in) :: alpha

!    LOG_DEBUG("multiplying matrix by scalar "//to_string(alpha))
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
    REAL(Kind(0d0)) :: TInitial, TTotal

    IF(.NOT.ASSOCIATED(qA)) then
!       LOG_DEBUG("qA not associated")
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


  !> Trace for QuTree: @f$ a = \mathrm{Tr} [ A ] @f$.
  !!
  !! @param qA Pointer to matrix A.
  !!
  !! @return The trace.
  FUNCTION SpAMM_Trace_QuTree (qA) RESULT(a)

    TYPE(QuTree), POINTER, INTENT(IN) :: qA
    REAL(SpAMM_KIND)                  :: a

    REAL(Kind(0d0)) :: TInitial, TTotal

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

    type(spamm_matrix_order_2), pointer, intent(in) :: A
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
    REAL(Kind(0d0)) :: TInitial, TTotal

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

    type(spamm_matrix_order_2), pointer, intent(in) :: A, B
    real(spamm_kind), intent(in), optional :: tolerance
    real(spamm_kind) :: trace_ab

    real(spamm_kind) :: local_tolerance

    trace_ab = 0

    if(.not. associated(A) .or. .not. associated(B)) then
!       LOG_INFO("either A or B are not associated")
       return
    endif

    if(.not. associated(A%root) .or. .not. associated(B%root)) then
!       LOG_INFO("either A%root or B%root are not associated")
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
    REAL(Kind(0d0))     :: TInitial, TTotal
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
    REAL(Kind(0d0)) :: TInitial, TTotal

    TInitial = SpAMM_Get_Time()
    !$OMP TASK SHARED(Norm,qA)
    Norm = SpAMM_Norm_Reduce_QuTree_Recur(qA)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_Get_Time()-TInitial
    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Norm_Reduce_QuTree",9)

  END FUNCTION SpAMM_Norm_Reduce_QuTree


  !> Inner product: \f$ (A,B) \f$.
  FUNCTION SpAMM_Dot_Product_BiTree(bA,bB) RESULT(dot)

    INTEGER              :: Depth
    TYPE(BiTree),POINTER :: bA,bB
    REAL(SpAMM_KIND)     :: dot
    REAL(Kind(0d0))                                  :: TInitial, TTotal
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
    REAL(Kind(0d0)) :: TInitial, TTotal

    IF(PRESENT(tolerance))THEN
       local_tolerance = tolerance
    ELSE
       local_tolerance = 0
    ENDIF

    TInitial = SpAMM_Get_Time()

    CALL SpAMM_Multiply_BiTree_x_Scalar(bC, 0d0)
    !$OMP TASK UNTIED SHARED(qA,bB,bC)
    CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC, qA, bB, local_tolerance)
    ! <<<<<<<<<<< This NORM should be done properly on the fly, in the recursive multiply >>>>>>>>>>>>>>>>>>>
    bC%Norm=NORM(bC)
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
    REAL(Kind(0d0)) :: TInitial, TTotal

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
    REAL(Kind(0d0))   :: TInitial, TTotal

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
#ifdef SPAMM_STORE_TRANSPOSE
    integer :: i, j
#endif

    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB)) THEN
!       LOG_DEBUG("qA: "//to_string(qA))
!       LOG_DEBUG("qB: "//to_string(qB))

       ! Apply the SpAMM condition.
       if(qA%Norm*qB%Norm <= threshold) then
!          LOG_DEBUG("going back up")
          RETURN
       endif

!#ifdef _OPENMP
!       CALL OMP_SET_LOCK(qC%lock)
!#endif
       IF(.NOT.ASSOCIATED(qC))THEN
          ! Allocate new node.
          CALL NewQuNode(qC, qA%i_lower, qB%j_lower, qA%i_upper, qB%j_upper)
       ENDIF
!#ifdef _OPENMP
!       CALL OMP_UNSET_LOCK(qC%lock)
!#endif

!       LOG_DEBUG("qC: "//to_string(qC))
!       LOG_DEBUG("   operations = "//to_string(qC%number_operations))

       ! At the bottom, calculate the product.
       IF(qC%i_upper-qC%i_lower+1 == SPAMM_BLOCK_SIZE) then
!#ifdef _OPENMP
!          CALL OMP_SET_LOCK(qC%lock)
!#endif
          IF(.NOT.ALLOCATED(qC%Blok))THEN
             ! Allocate new block.
             ALLOCATE(qC%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
             qC%Blok = 0
          ENDIF
!#ifdef _OPENMP
!          CALL OMP_UNSET_LOCK(qC%lock)
!#endif

!#if defined(_OPENMP) && ! defined(BIGLOCK)
!     CALL OMP_SET_LOCK(qC%lock)
!#endif
#ifdef SPAMM_STORE_TRANSPOSE
          do i = 1, SPAMM_BLOCK_SIZE
             do j = 1, SPAMM_BLOCK_SIZE
                qC%blok = qC%blok + dot_product(qA%blok(:, i), qB%blok(:, j))
             enddo
          enddo
#else
          qC%Blok = qC%Blok + matmul(qA%Blok, qB%Blok)
#endif
#ifdef SPAMM_COUNTERS
          qC%number_operations = qC%number_operations+SPAMM_BLOCK_SIZE**3
          qC%number_nonzeros = count_nonzero(qC%blok)
#endif

!#if defined(_OPENMP) && ! defined(BIGLOCK)
!     CALL OMP_UNSET_LOCK(qC%lock)
!#endif
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
!       LOG_DEBUG("either A or B is not associated")
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur

  !> Recursive part of scalar multiply with quadtree matrix, @f$ A \leftarrow \alpha A @f$.
  !!
  !! @note The number_operations counter is disabled.
  !!
  !! @param qA Pointer to quadtree.
  !! @param alpha The scalar.
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)

    TYPE(QuTree), POINTER :: qA
    REAL(SpAMM_KIND) :: alpha

    IF(.NOT.ASSOCIATED(qA)) then
!       LOG_DEBUG("qA not associated")
       RETURN
    endif

!    LOG_DEBUG("qA = "//to_string(qA))

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
       ! At the bottom, multiply the block.
       qA%Norm = qA%Norm*ABS(alpha)
       qA%Blok = alpha*qA%Blok
!#ifdef SPAMM_STORE_TRANSPOSE
!       qA%transpose_block = alpha*qA%transpose_block
!#endif
#ifdef SPAMM_COUNTERS
       ! Disable for now...
       !qA%number_operations = SPAMM_BLOCK_SIZE**2
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

    Trace = 0

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) THEN
       Trace = 0
       IF(.NOT.ASSOCIATED(qA)) RETURN
       DO I = 1, SPAMM_BLOCK_SIZE
          Trace = Trace+qA%Blok(I,I)
       ENDDO
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11).AND. &
         .NOT.ASSOCIATED(qA%Quad22))THEN
       Trace = 0
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

  !> Recursive part of trace for quadtree product, @f$ \mathrm{Tr} \left[ A
  !! \cdot B \right] @f$.
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

    Trace = 0

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
       Norm = 0
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

  !> Recursive linear algebra routines on row tree vectors: @f$ C \leftarrow A
  !! \dot B @f$.
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
!       LOG_DEBUG("qA: "//to_string(qA))
!       LOG_DEBUG("bB: "//to_string(bB))

       if(qA%Norm*bB%Norm < tolerance) then
!          LOG_DEBUG("norm_A*norm_B = "//to_string(qA%Norm*bB%Norm)//" < "//to_string(tolerance))
          return
       endif

       IF(.NOT.ASSOCIATED(bC))THEN
!          LOG_DEBUG("allocating new node in C bitree")
          call new(bC, bB%i_lower, bB%i_upper)
       ENDIF

       ! Blocks
       !      IF(bC%i_upper-bC%i_lower+1 == SPAMM_BLOCK_SIZE) THEN

       IF(ALLOCATED(qA%Blok).AND.ALLOCATED(bB%Vect))THEN
          ! Allocate
          IF(.NOT.ALLOCATED(bC%Vect))THEN
!             LOG_DEBUG("allocating new vec in C bitree")
             !$OMP CRITICAL
             ALLOCATE(bC%Vect(SPAMM_BLOCK_SIZE))
             bC%Vect = 0
             !$OMP END CRITICAL
          END IF
          ! Accumulate
          bC%Vect(1:SPAMM_BLOCK_SIZE)=bC%Vect+matmul(qA%Blok, bB%Vect)
       ELSE
!          LOG_DEBUG("descending")
          bC%Norm = 0
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
!          LOG_DEBUG("ascending")
       ENDIF
    else
!       LOG_DEBUG("either qA or bB are not associated")
    ENDIF
  END SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree_Recur


  !> \f$ L_2 \f$ norm for BiTrees.
  recursive function SpAMM_Norm_Reduce_BiTree_Recur(bA,Depth) result(Norm)

    type(BiTree), pointer :: bA
    integer               :: Depth
    real(SpAMM_KIND)      :: Norm, Norm0, Norm1

    if(.not.associated(bA))then
       Norm = 0
       return
    elseif(allocated(bA%Vect))then
       Norm=sum(bA%Vect(1:SPAMM_BLOCK_SIZE)**2)
       bA%Norm=sqrt(Norm)
    else
       Norm0 = 0
       Norm1 = 0

       !$OMP TASK UNTIED SHARED(bA) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Norm0=SpAMM_Norm_Reduce_BiTree_Recur(bA%sect1,Depth+1)
       !$OMP END TASK
       !$OMP TASK UNTIED SHARED(bA) &
       !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
       Norm1=SpAMM_Norm_Reduce_BiTree_Recur(bA%sect2,Depth+1)
       !$OMP END TASK
       !$OMP TASKWAIT
       Norm=Norm0+Norm1
       bA%Norm=sqrt(Norm)
    endif
  end function SpAMM_Norm_Reduce_BiTree_Recur

  !> Dot product for BiTrees
  RECURSIVE FUNCTION SpAMM_Dot_Product_BiTree_Recur(bA,bB,Depth) RESULT(Dot)
    TYPE(BiTree), POINTER :: bA,bB
    INTEGER               :: Depth
    REAL(SpAMM_KIND)      :: Dot, Dot0, Dot1

    Dot = 0

    IF(.NOT.ASSOCIATED(bA))THEN
       RETURN
    ELSEIF(.NOT.ASSOCIATED(bB))THEN
       RETURN
    ELSEIF(bA%Norm*bB%Norm<1D-20)THEN  !! <<<<<<<<<<<<  Needs to have a threshold passed in !! >>>>>>>>>>>>>
       RETURN
    ELSEIF(allocated(bA%Vect).AND.allocated(bB%Vect))THEN
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

!    LOG_DEBUG("q: "//to_string(bA%i_lower)//", "//to_string(bA%i_upper))

    IF(bA%i_upper-bA%i_lower+1 == SPAMM_BLOCK_SIZE .AND. ALLOCATED(bA%Vect))THEN
       bA%Norm=bA%Norm*ABS(alpha)
       bA%Vect=bA%Vect*alpha
!       LOG_DEBUG("vect = "//to_string(bA%vect(1)))
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

  !> Multiply two 2nd order matrices: \f$ C \leftarrow \alpha A B + \beta C
  !! \f$.
  !!
  !! If the tolerance is not given, then \f$ \tau = 0 \f$ is used, i.e. the
  !! product reverts to an exact dense product.
  !!
  !! @bug The code for \f$ \alpha \neq 1 \f$ and \f$ \beta \neq 0 \f$ is not
  !! implemented yet.
  !!
  !! @param A The matrix \f$ A \f$.
  !! @param B The matrix \f$ B \f$.
  !! @param C The matrix \f$ C \f$.
  !! @param tolerance The SpAMM tolerance \f$ \tau \f$.
  !! @param alpha The scalar \f$ \alpha \f$.
  !! @param beta The scalar \f$ \beta \f$.
  subroutine spamm_multiply_2nd_order_x_2nd_order (A, B, C, tolerance, alpha, beta)

    type(spamm_matrix_order_2), pointer, intent(in) :: A, B
    type(spamm_matrix_order_2), pointer, intent(inout) :: C
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
!       LOG_INFO("either A or B are not allocated")
       return
    endif

    if(associated(A, C)) then
!       LOG_FATAL("A == C")
       error stop
    endif

    if(associated(B, C)) then
!       LOG_FATAL("B == C")
       error stop
    endif

    if(.not. associated(C)) then
       call new(A%M, B%N, C)
    endif

!    LOG_DEBUG("multiplying A*B with tolerance "//to_string(local_tolerance))

    call reset_counters(C)
    call spamm_multiply_qutree_x_qutree(A%root, B%root, C%root, local_tolerance)

    if(associated(C%root)) then
       C%norm = C%root%norm
       C%number_operations = C%root%number_operations
       C%number_nonzeros = C%root%number_nonzeros
    endif

  end subroutine spamm_multiply_2nd_order_x_2nd_order

  subroutine spamm_multiply_2nd_order_x_1st_order (A, B, C, threshold)

    type(spamm_matrix_order_2), pointer, intent(in) :: A
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

!    LOG_DEBUG("multiply matrix*vector")
!    LOG_DEBUG("A: "//to_string(A))
!    LOG_DEBUG("B: "//to_string(B))
    call spamm_multiply_qutree_x_bitree(A%root, B%root, C%root, local_threshold)
  end subroutine spamm_multiply_2nd_order_x_1st_order

  !> Frobenius norm of 2nd order matrix. This function updates the norm on the
  !! matrix and returns the square of the norm.
  !!
  !! @param A The matrix.
  !!
  !! @return The squared Frobenius norm.
  function spamm_norm_reduce_matrix_2nd_order (A) result (norm)

    real(spamm_kind) :: norm
    type(spamm_matrix_order_2), pointer, intent(in) :: A

    A%norm = spamm_norm_reduce_qutree_recur(A%root)
    norm = A%norm

  end function spamm_norm_reduce_matrix_2nd_order
