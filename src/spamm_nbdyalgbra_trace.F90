module spamm_nbdyalgbra_trace

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals

  implicit none

CONTAINS

  FUNCTION SpAMM_tree_2d_symm_trace (a) RESULT(tr)

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN) :: A
    REAL(SpAMM_KIND)                              :: tr

    tr = SpAMM_tree_2d_symm_trace_recur(a)

  END FUNCTION SpAMM_tree_2d_symm_trace


  RECURSIVE FUNCTION SpAMM_tree_2d_symm_trace_recur(a) RESULT(tr)

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN) :: a
    REAL(SpAMM_KIND)                              :: tr
    INTEGER, DIMENSION(:,:),  pointer             :: bb
    integer                                       :: i

    tr=SpAMM_Zero
    IF(.NOT.ASSOCIATED(a))RETURN

    bb => a%frill%bndbx 

    IF( bb(1,1)-bb(0,1) == SBS )THEN ! Leaf condition ? 

       do i = 1,SBS
          tr=tr+a%chunk(i,i)
       enddo

    ELSE

       tr=tr+SpAMM_tree_2d_symm_trace_recur(a%child_00)
       tr=tr+SpAMM_tree_2d_symm_trace_recur(a%child_01)
       tr=tr+SpAMM_tree_2d_symm_trace_recur(a%child_11)

    ENDIF

  END FUNCTION SpAMM_tree_2d_symm_trace_recur
















!!$
!!$
!!$  FUNCTION SpAMM_Trace_QuTree_Product(qA, qB, threshold) RESULT(a)
!!$
!!$    TYPE(QuTree), POINTER      :: qA,qB
!!$    REAL(SpAMM_KIND), OPTIONAL :: threshold
!!$    REAL(SpAMM_KIND)           :: a
!!$
!!$    INTEGER            :: Depth
!!$    REAL(SpAMM_KIND)   :: multiplyThreshold
!!$    REAL(Kind(0d0)) :: TInitial, TTotal
!!$
!!$    Depth=0
!!$
!!$    IF(PRESENT(threshold)) THEN
!!$       multiplyThreshold = threshold
!!$    ELSE
!!$       multiplyThreshold = 0
!!$    ENDIF
!!$
!!$    TInitial = SpAMM_Get_Time()
!!$
!!$    !$OMP TASK UNTIED SHARED(qA,a)
!!$    a = SpAMM_Trace_QuTree_Product_Recur(qA, qB, multiplyThreshold)
!!$    !$OMP END TASK
!!$    !$OMP TASKWAIT
!!$
!!$    TTotal=SpAMM_Get_Time()-TInitial
!!$    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree_Product",7)
!!$
!!$  END FUNCTION SpAMM_Trace_QuTree_Product
!!$
!!$  !> Trace for matrix product, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
!!$  !!
!!$  !! @param A Pointer to matrix A.
!!$  !! @param B Pointer to matrix B.
!!$  !!
!!$  !! @return The trace of the matrix produce, @f$ \mathrm{Tr} \left[ A \cdot B \right] @f$.
!!$  function spamm_trace_2nd_order_product (A, B, tolerance) result (trace_ab)
!!$
!!$    type(spamm_matrix_order_2), pointer, intent(in) :: A, B
!!$    real(spamm_kind), intent(in), optional :: tolerance
!!$    real(spamm_kind) :: trace_ab
!!$
!!$    real(spamm_kind) :: local_tolerance
!!$
!!$    trace_ab = 0
!!$
!!$    if(.not. associated(A) .or. .not. associated(B)) then
!!$!       LOG_INFO("either A or B are not associated")
!!$       return
!!$    endif
!!$
!!$    if(.not. associated(A%root) .or. .not. associated(B%root)) then
!!$!       LOG_INFO("either A%root or B%root are not associated")
!!$       return
!!$    endif
!!$
!!$    if(present(tolerance)) then
!!$       local_tolerance = tolerance
!!$    else
!!$       local_tolerance = 0
!!$    endif
!!$
!!$    trace_ab = spamm_trace_qutree_product(A%root, B%root, local_tolerance)
!!$
!!$  end function spamm_trace_2nd_order_product
!!$
!!$


end module spamm_nbdyalgbra_trace


