module spamm_nbdyalgbra_plus

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals
  use spamm_nbdyalgbra_times 
  implicit none

CONTAINS

  !++NBODYADDTN: SpAMM generalized n-body algebras for the add ____ NBODYADDTN _________________
  !++NBODYADDTN: generalized addition (+)
  !++NBODYADDTN:   ... [TREE-ONE-D + TREE-ONE-D] ... [TREE-ONE-D + TREE-ONE-D] ...   
  !++NBODYADDTN:   SpAMM_tree_1d_plus_tree_1d 


  !++NBODYADDTN:     c => inplace*c + beta*b
  FUNCTION SpAMM_tree_1d_plus_tree_1d (C, beta, B, inplace) RESULT(D)   

    TYPE(SpAMM_tree_1d), POINTER, INTENT(IN)    :: B
    TYPE(SpAMM_tree_1d), POINTER, INTENT(INOUT) :: C
    REAL(SpAMM_KIND),             INTENT(IN)    :: Beta
    real(spamm_kind),    OPTIONAL,intent(in)    :: inplace
    TYPE(SpAMM_tree_1d), POINTER                :: d
    REAL(SpAMM_KIND)                            :: Alpha
    !
    D=>C
    !
    if(.not. associated(B))RETURN
    if(beta==SpAMM_zero)   RETURN
    !
    IF(PRESENT(inplace))THEN
       alpha=inplace
    ELSE
       alpha=SpAMM_one
    ENDIF

    CALL SpAMM_tree_1d_plus_tree_1d_inplace_recur(d, B, alpha, beta )
    
  END FUNCTION SpAMM_tree_1d_plus_tree_1d 

  ! for tree_1d, A = alpha*A + beta*B
  RECURSIVE SUBROUTINE SpAMM_tree_1d_plus_tree_1d_inplace_recur(a, b, alpha, beta)

    TYPE(SpAMM_tree_1d), POINTER                :: A
    TYPE(SpAMM_tree_1d), POINTER, INTENT(IN)    :: B
    REAL(SpAMM_KIND),             INTENT(IN)    :: alpha, beta
    logical                                     :: TA, TB

    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)

    IF(.NOT.TA.AND..NOT.TB)THEN

       RETURN

    ELSEIF(TA .AND. .NOT.TB) THEN

       ! a=alpha*a
       CALL SpAMM_scalar_times_tree_1d_recur(alpha, a)

    ELSEIF(TA.AND.TB)THEN

       IF(b%frill%leaf)then
          ! A = alpha*A + beta*B

          a%frill%init=.FALSE.
          a%chunk(1:SBS) = alpha*a%chunk(1:SBS) + beta*b%chunk(1:SBS)
          a%frill%flops = a%frill%flops + 3*SBS                  

       ELSE ! recursively descend ...

          CALL SpAMM_tree_1d_plus_tree_1d_inplace_recur(SpAMM_construct_tree_1d_0(a), & !0>
                                                        b%child_0, alpha, beta)

          CALL SpAMM_tree_1d_plus_tree_1d_inplace_recur(SpAMM_construct_tree_1d_1(a), & !1> 
                                                        b%child_1, alpha, beta)
       ENDIF

    ENDIF

    CALL SpAMM_redecorate_tree_1d(a)

  END SUBROUTINE SpAMM_tree_1d_plus_tree_1d_inplace_recur


  !!
  !! ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ...
  !!


  FUNCTION SpAMM_scalar_plus_tree_2d_symm(alpha, a) RESULT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER  :: a 
    TYPE(SpAMM_tree_2d_symm), POINTER  :: d 
    REAL(SpAMM_KIND),  INTENT(IN) :: alpha

    if(.not.associated(a))then
       d=>null()
       return
    else
       d=>a
       call SpAMM_scalar_plus_tree_2d_symm_recur(alpha, d)
    endif

  END FUNCTION SpAMM_scalar_plus_tree_2d_symm

  recursive subroutine SpAMM_scalar_plus_tree_2d_symm_recur(alpha, a)

   TYPE(SpAMM_tree_2d_symm), POINTER  :: a 
   REAL(SpAMM_KIND),  INTENT(IN) :: alpha
   INTEGER, dimension(1:2)       :: lo,hi
   INTEGER                       :: i

   if(.not.associated(a))return

   IF(a%frill%leaf)THEN

      a%frill%init=.FALSE.
      lo=a%frill%bndbx(0,:)
      hi=a%frill%bndbx(1,:)
      
      do i=1,hi(1)-lo(1)+1
         a%chunk(i,i)=a%chunk(i,i)+alpha
      enddo

    ELSE
       ! child along [00]:
       CALL SpAMM_scalar_plus_tree_2d_symm_recur(alpha, SpAMM_construct_tree_2d_symm_00(a))
       ! child along [11]:
       CALL SpAMM_scalar_plus_tree_2d_symm_recur(alpha, SpAMM_construct_tree_2d_symm_11(a))
    ENDIF
    
    CALL SpAMM_redecorate_tree_2d_symm(a)

  end subroutine SpAMM_scalar_plus_tree_2d_symm_recur


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

    IF(PRESENT(C))THEN ! we are going for an in place plus with an existing C:

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

       ! a = b (zero thresholding)
       CALL SpAMM_tree_2d_symm_copy_tree_2d_symm_recur(a, b, SpAMM_zero )

       ! a = beta*b = beta*a  
!       CALL SpAMM_scalar_times_tree_2d_symm_recur(beta, a)

    ELSEIF(TA .AND. .NOT.TB) THEN

       ! a=alpha*a
!       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a)

    ELSEIF(TA.AND.TB)THEN

       IF(b%frill%leaf)then

          ! A = alpha*A + beta*B
          a%frill%init=.false.
          a%chunk(1:SBS,1:SBS) = alpha*a%chunk(1:SBS,1:SBS) + beta*b%chunk(1:SBS,1:SBS)
          a%frill%flops = a%frill%flops + 3*SBS2                  

       ELSE

          ! recursively decend, possibly building out A if nessesary ...
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_00(a), & !00>
                                                                  b%child_00, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_01(a), & !01> 
                                                                  b%child_01, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_10(a), & !10> 
                                                                  b%child_10, alpha, beta)
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
          c%frill%init=.FALSE.
          c%chunk(1:sbs,1:sbs)=c%chunk(1:sbs,1:sbs)+alpha*a%chunk(1:sbs,1:sbs)+beta*b%chunk(1:sbs,1:sbs)
          c%frill%flops=c%frill%flops+3*sbs2

       ELSE
          ! recursively decend, popping new children as needed ...
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_00(c), & !00> 
                                                           a%child_00, b%child_00, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_01(c), & !01>
                                                           a%child_01, b%child_01, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_10(c), & !10>
                                                           a%child_10, b%child_10, alpha, beta)
          CALL SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_11(c), & !11>
                                                           a%child_11, b%child_11, alpha, beta)
       ENDIF

       CALL SpAMM_redecorate_tree_2d_symm(c)

    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_plus_tree_2d_symm_recur

END module spamm_nbdyalgbra_plus
