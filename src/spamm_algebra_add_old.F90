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

    type(spamm_matrix_order_2), pointer, intent(inout) :: A
    type(spamm_matrix_order_2), pointer, intent(in) :: B
    real(spamm_kind), intent(in), optional :: alpha, beta

!    LOG_DEBUG("Adding matrices: alpha = "//to_string(alpha)//", beta = "//to_string(beta))

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
    REAL(Kind(0d0)) :: TInitial, TTotal

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

!    LOG_DEBUG("Adding "//to_string(inplace_alpha)//"*A + "//to_string(inplace_beta)//"*B")

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

    type(spamm_matrix_order_2), pointer, intent(inout) :: A
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
    REAL(Kind(0d0)) :: TInitial, TTotal

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


  !> Add to binary trees: \f$ C \leftarrow \alpha A + \beta B \f$.
  SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace(bC,Alpha,bA,Beta,bB)

    TYPE(BiTree), POINTER, INTENT(IN)    :: bA,bB
    TYPE(BiTree), POINTER, INTENT(INOUT) :: bC
    REAL(SpAMM_KIND), OPTIONAL           :: Alpha,Beta
    INTEGER                              :: Depth
    REAL(Kind(0d0))                   :: TInitial, TTotal

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
    bC%Norm=NORM(bC)
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
    REAL(Kind(0d0))                   :: TInitial, TTotal

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
    bA%Norm=NORM(bA)
  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace


  !> Add 2 quadtree matrices, @f$ A \leftarrow \alpha A + \beta B @f$.
  !!
  !! @param qA [inout] Pointer to matrix A.
  !! @param qB [in] Pointer to matrix B.
  !! @param alpha Factor @f$ \alpha @f$.
  !! @param beta Factor @f$ \beta @f$.
  RECURSIVE SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA, qB, alpha, beta)

    type(QuTree), pointer :: qA,qB
    real(SpAMM_KIND) :: alpha, beta
    logical :: TA, TB

    TA=ASSOCIATED(qA)
    TB=ASSOCIATED(qB)

    IF(TA.AND.TB)THEN
!       LOG_DEBUG("qA = "//to_string(qA))
!       LOG_DEBUG("qB = "//to_string(qB))

       IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) THEN
          qA%Blok = alpha*qA%Blok+beta*qB%Blok
#ifdef SPAMM_STORE_TRANSPOSE
!          qA%transpose_block = alpha*qA%transpose_block+beta*qB%transpose_block
#endif
          qA%norm = sqrt(sum(qA%blok**2))
#ifdef SPAMM_COUNTERS
          qA%number_nonzeros = count_nonzero(qA%blok)
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
!       LOG_DEBUG("either A or B are not associated")
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

!    LOG_DEBUG("q:"//to_string(qA%i_lower)//" "//to_string(qA%j_lower))
!    LOG_DEBUG("  "//to_string(qA%i_upper)//" "//to_string(qA%j_upper))

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
       if(.not. allocated(qA%blok)) then
          ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
          qA%Blok = 0
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

  !> BiTree in place add: \f$ A \leftarrow \alpha A + \beta B \f$.
  RECURSIVE SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA,bB,Depth)
    TYPE(BiTree),POINTER :: bA,bB
    INTEGER              :: Depth
    LOGICAL              :: TA, TB


    TA=ASSOCIATED(bA)
    TB=ASSOCIATED(bB)
    IF(TA.AND.TB)THEN
       IF(ALLOCATED(bA%Vect).AND.ALLOCATED(bB%Vect))THEN
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
       IF(ALLOCATED(bA%Vect).AND.ALLOCATED(bB%Vect))THEN
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
