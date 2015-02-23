module spamm_nbdyalgbra_plus

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals

  implicit none

CONTAINS

  FUNCTION SpAMM_tree_2d_symm_plus_tree_2d_symm (A, B, alpha, beta, C) RESULT(D)   

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(INOUT)           :: A, B
    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(INOUT), OPTIONAL :: C
    real(spamm_kind), intent(in), optional                     :: alpha, beta

    TYPE(SpAMM_tree_2d_symm), POINTER                          :: d

    REAL(SpAMM_KIND)                                           :: Local_Alpha,Local_Beta
    !
    D=>NULL()
    !
    if(.not. associated(A))RETURN
    if(.not. associated(B))RETURN
    !
    IF(PRESENT(Alpha))THEN; Local_Alpha=Alpha; ELSE; Local_Alpha=SpAMM_One; ENDIF
    IF(PRESENT(Beta ))THEN; Local_Beta =Beta;  ELSE; Local_Beta=SpAMM_One;  ENDIF

    IF(PRESENT(C))THEN ! we are going for an in place add with an existing C:

       IF(ASSOCIATED(B,C))THEN  ! if passed in C is B, then in place accumulation on B ...

          ! B => alpha*A + beta*B 
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(B, A, Local_beta, Local_alpha)
          D=>B

       ELSEIF(ASSOCIATED(A,C))THEN ! if passed in C is A, then in place accumulation on A ...

          ! A => alpha*A + beta*B ...
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(A, B, Local_alpha, Local_beta)
          D=>A

       ELSE  ! C is passed in as a seperate channel for accumulation ... 

          ! C => C + alpha*A + beta*B
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(C, A, B, Local_alpha, Local_beta)
          D=>C

       ENDIF

    ELSE  ! We need a clean tree_2d_symm at this point ...

       D => SpAMM_new_top_tree_2d_symm(A%frill%NDimn)
       ! D => D + alpha*A + beta*B
       CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(C, A, B, Local_alpha, Local_beta)
       D=>C

    ENDIF

     CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(D, A, B, alpha, beta)

  END FUNCTION SpAMM_tree_2d_symm_plus_tree_2d_symm 

  ! for tree_2d_symm, A = A + alpha*A + beta*B

  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(a, b, alpha, beta)

    TYPE(SpAMM_tree_2d_symm), POINTER                :: A
    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: B
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    logical                                          :: TA, TB

    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)

    IF(.NOT.TA.AND.TB)THEN

       ! a = b
       CALL SpAMM_tree_2d_symm_copy_tree_2d_symm_recur(a, b)

       ! a = beta*b = beta*a  
       CALL SpAMM_scalar_times_tree_2d_symm_recur(beta, a)

    ELSEIF(TA .AND. .NOT.TB) THEN

       ! a=alpha*a
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a)

    ELSEIF(TA.AND.TB)THEN

       IF(b%frill%leaf)then

          ! A = alpha*A + beta*B
          a%chunk(1:SBS,1:SBS) = alpha*a%chunk(1:SBS,1:SBS) + beta*b%chunk(1:SBS,1:SBS)
          a%frill%flops = a%frill%flops + 3*SBS2                  

       ELSE

          ! recursively decend, possibly building out A if nessesary ...
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_00(a), & !00>
                                                                  b%child_00, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_01(a), & !01> 
                                                                  b%child_01, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_11(a), & !11> 
                                                                  b%child_11, alpha, beta)
       ENDIF

       ! enrich the resultnt
       CALL SpAMM_redecorate_tree_2d_symm(a)

    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur

  ! for tree_2d_symm: C = C + alpha*A+beta*B
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(C, A, B, alpha, beta)

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A,B
    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
    REAL(SpAMM_KIND)                                 :: alpha, beta
    logical                                          :: TA, TB

    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)

    IF(TA .AND. .NOT.TB) THEN

       ! C = C + alpha*A
       CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(c, a, alpha, SpAMM_One)

    ELSEIF(.NOT.TA.AND.TB)THEN

       ! C = C + beta*B
       CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(c, b, SpAMM_One, beta)

    ELSEIF(TA.AND.TB)THEN

       ! leaf situation ...
       IF(c%frill%leaf)THEN

          ! c = c + alpha*a + beta*b
          c%chunk(1:sbs,1:sbs)=c%chunk(1:sbs,1:sbs)+alpha*a%chunk(1:sbs,1:sbs)+beta*b%chunk(1:sbs,1:sbs)
          c%frill%flops=c%frill%flops+3*sbs2

       ELSE

          ! recursively decend, popping new children as needed ...
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_00(c), & !00> 
                                                           a%child_00, b%child_00, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_01(c), & !01>
                                                           a%child_01, b%child_01, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_11(c), & !11>
                                                           a%child_11, b%child_11, alpha, beta)
       ENDIF

       CALL SpAMM_redecorate_tree_2d_symm(c)

    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_plus_tree_2d_symm_recur

!!$  FUNCTION SpAMM_Add_Tree1d_2_Tree1d (A, B, alpha, beta, C) RESULT(D)   
!!$    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT)           :: A, B
!!$    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT), OPTIONAL :: C
!!$    TYPE(SpAMM_Tree_1d), POINTER                          :: D
!!$    real(spamm_kind), intent(in), optional                :: alpha, beta
!!$    REAL(SpAMM_KIND)                                      :: Local_Alpha,Local_Beta
!!$    !
!!$    D=>NULL()
!!$    !
!!$    ! hmmm... this needs some thought ...
!!$    !if(.not. associated(A))RETURN
!!$    !if(.not. associated(B))RETURN
!!$    !
!!$    IF(PRESENT(Alpha))THEN
!!$       Local_Alpha=Alpha
!!$    ELSE
!!$       Local_Alpha=SpAMM_One
!!$    ENDIF
!!$    IF(PRESENT(Beta))THEN
!!$       Local_Beta=Beta
!!$    ELSE
!!$       Local_Beta=SpAMM_One
!!$    ENDIF
!!$
!!$    IF(PRESENT(C))THEN
!!$       !  C=>Add(A,B,a,C): in place on C, C can be A, B or other.
!!$       IF(ASSOCIATED(B,C))THEN
!!$          ! B=>a*A+b*B
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(B, A, Local_beta, Local_alpha)
!!$          D=>B
!!$       ELSEIF(ASSOCIATED(A,C))THEN
!!$          ! A=>a*A+b*B
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A, B, Local_alpha, Local_beta)
!!$          D=>A
!!$       ELSE
!!$          ! C=>a*A+B
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(C, A, B, Local_alpha, Local_beta)
!!$          D=>C
!!$       ENDIF
!!$    ELSE
!!$       ! D=>a*A+B
!!$       !       CALL New(D, A%decoration )
!!$       CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(D, A, B, Local_alpha, Local_beta)       
!!$    ENDIF
!!$  END FUNCTION SpAMM_Add_Tree1d_2_Tree1d
!!$
!!$  RECURSIVE SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A, B, alpha, beta)
!!$    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT) :: A
!!$    TYPE(SpAMM_Tree_1d), POINTER, INTENT(IN)    :: B
!!$    REAL(SpAMM_KIND),             INTENT(IN)    :: alpha, beta
!!$    logical                                     :: TA, TB
!!$    TA=ASSOCIATED(A)
!!$    TB=ASSOCIATED(B)
!!$    IF(.NOT.TA.AND.TB)THEN
!!$       CALL SpAMM_Copy_Tree1d_2_Tree1d(B, A)
!!$       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, beta)
!!$    ELSEIF(TA .AND. .NOT.TB) THEN
!!$       ! Multiply A tree with alpha.
!!$       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, alpha)
!!$    ELSEIF(TA.AND.TB)THEN
!!$       IF(ALLOCATED(A%chunk).AND.ALLOCATED(B%chunk))THEN
!!$          A%chunk=alpha*A%chunk+beta*B%chunk
!!$       ELSE
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A%child_0, B%child_0, alpha, beta)
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A%child_1, B%child_1, alpha, beta)
!!$       ENDIF
!!$       CALL ReDecorate(A)
!!$    ENDIF
!!$  END SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur
!!$
!!$  RECURSIVE SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_Recur(A, B, C, alpha, beta)
!!$    TYPE(SpAMM_Tree_1d), POINTER, INTENT(IN)    :: A,B
!!$    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT) :: C
!!$    REAL(SpAMM_KIND)                            :: alpha, beta
!!$    logical                                     :: TA, TB
!!$    TA=ASSOCIATED(A)
!!$    TB=ASSOCIATED(B)
!!$    IF(.NOT.TA.AND.TB)THEN
!!$       CALL SpAMM_Copy_Tree1d_2_Tree1d(B, A)
!!$       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, beta)
!!$    ELSEIF(TA .AND. .NOT.TB) THEN
!!$       ! Multiply A tree with alpha.
!!$       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, alpha)
!!$    ELSEIF(TA.AND.TB)THEN
!!$       IF(ALLOCATED(A%chunk).AND.ALLOCATED(B%chunk))THEN
!!$          A%chunk=alpha*A%chunk+beta*B%chunk
!!$       ELSE
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(A%child_0, B%child_0, C%child_0, alpha, beta)
!!$          CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(A%child_1, B%child_1, C%child_1, alpha, beta)
!!$       ENDIF
!!$       CALL ReDecorate(C)
!!$    ENDIF
!!$  END SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_Recur
!!$

!  INCLUDE 'spamm_algebra_add_old.F90'

END module spamm_nbdyalgbra_plus
