
  FUNCTION SpAMM_Add_Tree_2d_2_Tree_2d (A, B, alpha, beta, C) RESULT(D)   
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT)           :: A, B
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT), OPTIONAL :: C
    TYPE(SpAMM_Tree_2d), POINTER                          :: D=>NULL()
    real(spamm_kind), intent(in), optional                :: alpha, beta
    if(.not. associated(A))RETURN
    if(.not. associated(B))RETURN
    IF(PRESENT(C))THEN
       !  C=>Add(A,B,a,C): in place on C, C can be A, B or other.
       IF(ASSOCIATED(B,C))THEN
          ! B=>a*A+b*B
          CALL SpAMM_Add_Tree_2d_2_Tree_2d_InPlace_Recur(B, A, beta, alpha)
          D=>B
       ELSEIF(ASSOCIATED(A,C))THEN
          ! A=>a*A+b*B
          CALL SpAMM_Add_Tree_2d_2_Tree_2d_InPlace_Recur(A, B, alpha, beta)
          D=>A
       ELSE
          ! C=>a*A+B
          CALL SpAMM_Add_Tree_2d_2_Tree_2d_Recur(C, A, B, alpha, beta)
          D=>C
       ENDIF
    ELSE
       ! D=>a*A+B
       CALL New(D, A%decoration )
       CALL SpAMM_Add_Tree_2d_2_Tree_2d_Recur(D, A, B, alpha, beta)       
    ENDIF
  END FUNCTION SpAMM_Add_Tree_2d_2_Tree_2d

  RECURSIVE SUBROUTINE SpAMM_Add_Tree_2d_2_Tree_2d_InPlace_Recur(A, B, alpha, beta)
    TYPE(SpAMM_Tree_2d), POINTER, pointer :: qA,qB
    REAL(SpAMM_KIND)                      :: alpha, beta
    logical :: TA, TB
    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)
    IF(.NOT.TA.AND.TB)THEN
       CALL SpAMM_Copy_Tree_2d_2_Tree_2d(B, A)
       CALL SpAMM_Multiply_Tree_2d_x_Scalar_recur(A, beta)
    ELSEIF(TA .AND. .NOT.TB) THEN
       ! Multiply A tree with alpha.
       CALL SpAMM_Multiply_Tree_2d_x_Scalar_recur(A, alpha)
    ELSEIF(TA.AND.TB)THEN
       IF(ALLOCATED(A%Blok).AND.ALLOCATED(B%Blok))THEN
          A%Blok=alpha*A%Blok+beta*B%Blok
       ELSE
          CALL SpAMM_Add_Tree_2d_2_Tree_2d_InPlace_Recur(A%child_00, B%child_00, alpha, beta)
          CALL SpAMM_Add_Tree_2d_2_Tree_2d_InPlace_Recur(A%child_01, B%child_01, alpha, beta)
          CALL SpAMM_Add_Tree_2d_2_Tree_2d_InPlace_Recur(A%child_11, B%child_11, alpha, beta)
       ENDIF
       CALL ReDecorate(A)
    ENDIF
  END SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur
