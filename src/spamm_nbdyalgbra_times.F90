module spamm_nbdyalgbra_times

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration

  implicit none

CONTAINS

  ! c => alpha*c + beta*(a.b)
  FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau2, alpha, beta, c) RESUT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha, beta
    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: C
    REAL(SpAMM_KIND)                                           :: Tau2
    TYPE(SpAMM_tree_2d_symm), POINTER                          :: D
    INTEGER                                                    :: Depth

    ! figure the starting conditions ...
    if(present(c))then
       d => c
    else
       d => NULL()
       if(present(alpha))write(*,*) ' multiplying though by ',alpha, ' whilst c uninited, dumbass ...'
    endif

    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! need a new tree, lets instantiate one ... 
    if(.not.associated(d)) &
       d => SpAMM_new_top_tree_2d_symm ( a%frill%ndimn(1), b%frill%ndimn(2) )

    if(present(alpha))then
       CALL SpAMM_scalar_times_tree_2d_symm(alpha, d)
    endif

    if(present(beta))then

    else
       CALL
    endif




!!$!    LOG_DEBUG("resetting C")
!!$    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, 0d0)
!!$
!!$!    LOG_DEBUG("recursive multiplication")
!!$
!!$    !$OMP TASK UNTIED SHARED(qA,qB,qC)
!!$    CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC, qA, qB, local_threshold)
!!$    !$OMP END TASK
!!$
!!$    !$OMP END MASTER
!!$
!!$    !$OMP END PARALLEL
!!$
!!$    qC%norm = norm(qC)
!!$    qC%norm = sqrt(qC%norm)
!!$
!!$    TTotal=SpAMM_Get_Time()-TInitial
!!$    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)
!!$
!!$  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree



  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur(C, A, B, Tau2, Depth )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(INOUT) :: C
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c11
    INTEGER, DIMENSION(:,:),  pointer                :: bb

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! spamm criteria
    if( a%Norm2 * b%Norm2 <= Tau2 ) return

    bb => c%frill%bndbx 
    IF( bb(1,1)-bb(0,1) == SBS )THEN ! Leaf condition ? 

       if(.not.allocated(c%chunk))then
          allocate(c%chunk(1:SBS,1:SBS))
          c%chunk=SpAMM_Zero
          c%frill%flops=0
       endif       

       c%chunk = c%chunk + matmul( a%chunk(1:SBS,1:SBS), b%chunk(1:SBS,1:SBS) )
       c%frill%flops = c%frill%flops + SBS3

    ELSE

       ! find some memory, new or old ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first tranch of 00,01 and 11 memory 
       CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(c00, a_child_00, b%child_00, Tau2, Depth+1)      
       CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(c01, a_child_00, b%child_01, Tau2, Depth+1)       
       CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(c11, a_child_10, b%child_01, Tau2, Depth+1)
       ! another pass with the same tranch
       CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(c00, a_child_01, b%child_10, Tau2, Depth+1)
       CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(c01, a_child_01, b%child_11, Tau2, Depth+1)
       CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(c11, a_child_11, b%child_11, Tau2, Depth+1)
       !
    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur


!!$  SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree(qA, qB, qC, threshold)
!!$
!!$    TYPE(QuTree), POINTER, INTENT(IN) :: qA, qB
!!$    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
!!$    REAL(SpAMM_KIND), OPTIONAL :: threshold
!!$
!!$    real(spamm_kind) :: local_threshold
!!$    REAL(kind(0d0)) :: TInitial, TTotal
!!$
!!$    TInitial = SpAMM_Get_Time()
!!$
!!$    if(.not. associated(qA) .or. .not. associated(qB)) then
!!$!       LOG_DEBUG("either A or B are not associated")
!!$       return
!!$    endif
!!$
!!$    if(.not. associated(qC)) then
!!$       CALL NewQuNode(qC, qA%i_lower, qB%j_lower, qA%i_upper, qB%j_upper)
!!$    endif
!!$
!!$    !$OMP PARALLEL
!!$
!!$    ! The master thread will lead execution of the product. All subsequent tasks
!!$    ! are untied and can be executed by any thread in the thread group.
!!$    !$OMP MASTER
!!$
!!$#ifdef _OPENMP
!!$!    LOG_INFO("Multiply on "//to_string(omp_get_num_threads())//" OpenMP threads")
!!$#endif
!!$
!!$    IF(PRESENT(threshold))THEN
!!$       local_threshold = threshold
!!$    ELSE
!!$       local_threshold = 0
!!$    ENDIF
!!$
!!$    qC%number_operations = 0
!!$
!!$!    LOG_DEBUG("resetting C")
!!$    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, 0d0)
!!$
!!$!    LOG_DEBUG("recursive multiplication")
!!$
!!$    !$OMP TASK UNTIED SHARED(qA,qB,qC)
!!$    CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC, qA, qB, local_threshold)
!!$    !$OMP END TASK
!!$
!!$    !$OMP END MASTER
!!$
!!$    !$OMP END PARALLEL
!!$
!!$    qC%norm = norm(qC)
!!$    qC%norm = sqrt(qC%norm)
!!$
!!$    TTotal=SpAMM_Get_Time()-TInitial
!!$    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)
!!$
!!$  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree
!!$
!!$  SUBROUTINE SpAMM_Multiply_QuTreeT_x_QuTree(qA, qB, qC, threshold)
!!$
!!$    TYPE(QuTree), POINTER, INTENT(IN) :: qA, qB
!!$    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
!!$    REAL(SpAMM_KIND), OPTIONAL :: threshold
!!$
!!$    real(spamm_kind) :: local_threshold
!!$    real(kind(0d0)) :: TInitial, TTotal
!!$
!!$    TInitial = SpAMM_Get_Time()
!!$
!!$    if(.not. associated(qA) .or. .not. associated(qB)) then
!!$!       LOG_DEBUG("either A or B are not associated")
!!$       return
!!$    endif
!!$
!!$    if(.not. associated(qC)) then
!!$       CALL NewQuNode(qC, qA%i_lower, qB%j_lower, qA%i_upper, qB%j_upper)
!!$    endif
!!$
!!$    !$OMP PARALLEL
!!$
!!$    ! The master thread will lead execution of the product. All subsequent tasks
!!$    ! are untied and can be executed by any thread in the thread group.
!!$    !$OMP MASTER
!!$
!!$#ifdef _OPENMP
!!$!    LOG_INFO("Multiply on "//to_string(omp_get_num_threads())//" OpenMP threads")
!!$#endif
!!$
!!$    IF(PRESENT(threshold))THEN
!!$       local_threshold = threshold
!!$    ELSE
!!$       local_threshold = 0
!!$    ENDIF
!!$
!!$    qC%number_operations = 0
!!$
!!$!    LOG_DEBUG("resetting C")
!!$    CALL SpAMM_Multiply_QuTree_x_Scalar(qC, 0d0)
!!$
!!$!    LOG_DEBUG("recursive multiplication")
!!$
!!$    !$OMP TASK UNTIED SHARED(qA,qB,qC)
!!$    !    CALL SpAMM_Multiply_QuTreeT_x_QuTree_Recur(qC, qA, qB, local_threshold)
!!$    !$OMP END TASK
!!$
!!$    !$OMP END MASTER
!!$
!!$    !$OMP END PARALLEL
!!$
!!$    qC%norm = norm(qC)
!!$    qC%norm = sqrt(qC%norm)
!!$
!!$    TTotal=SpAMM_Get_Time()-TInitial
!!$    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)
!!$
!!$  END SUBROUTINE SpAMM_Multiply_QuTreeT_x_QuTree
!!$
!!$  !> Scalar multiply of 2nd order matrix: @f$ A \leftarrow \alpha A @f$.
!!$  !!
!!$  !! @param A The matrix.
!!$  !! @param alpha The scalar \f$ \alpha \f$.
!!$  subroutine spamm_multiply_2nd_order_x_scalar (A, alpha)
!!$
!!$    type(spamm_matrix_order_2), pointer, intent(inout) :: A
!!$    real(spamm_kind), intent(in) :: alpha
!!$
!!$!    LOG_DEBUG("multiplying matrix by scalar "//to_string(alpha))
!!$    call spamm_multiply_qutree_x_scalar(A%root, alpha)
!!$    A%number_operations = A%root%number_operations
!!$
!!$  end subroutine spamm_multiply_2nd_order_x_scalar
!!$
!!$  !> Scalar multiply: @f$ A \leftarrow alpha A @f$.
!!$  !!
!!$  !! @param qA Pointer to matrix A.
!!$  !! @param alpha Scalar @f$ \alpha @f$.
!!$  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar(qA, alpha)
!!$
!!$    TYPE(QuTree), POINTER        :: qA
!!$    REAL(SpAMM_KIND), INTENT(IN) :: alpha
!!$
!!$    INTEGER            :: Depth
!!$    REAL(Kind(0d0)) :: TInitial, TTotal
!!$
!!$    IF(.NOT.ASSOCIATED(qA)) then
!!$!       LOG_DEBUG("qA not associated")
!!$       RETURN
!!$    endif
!!$
!!$    Depth=0
!!$    TInitial = SpAMM_Get_Time()
!!$
!!$    qA%number_operations = 0
!!$
!!$    !$OMP TASK UNTIED SHARED(qA)
!!$    CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA, alpha)
!!$    !$OMP END TASK
!!$
!!$    !$OMP TASKWAIT
!!$
!!$    qA%norm = norm(qA)
!!$    qA%norm = sqrt(qA%norm)
!!$
!!$    TTotal=SpAMM_Get_Time()-TInitial
!!$    !CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_Scalar",3)
!!$
!!$  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar




end module spamm_nbdyalgbra_times
