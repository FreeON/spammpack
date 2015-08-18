module spamm_algebra_add
  use spamm_types
  use spamm_globals
  use spamm_decorate
  use spamm_management
  use spamm_utilities
  implicit none
  PRIVATE
  PUBLIC :: Add
  !> Interface for additions operations between different SpAMM types.
  INTERFACE Add
!     MODULE PROCEDURE SpAMM_Add_Tree2dfull_2_Tree2dfull
     MODULE PROCEDURE SpAMM_Add_Tree2d_2_Tree2d
     MODULE PROCEDURE SpAMM_Add_Tree1d_2_Tree1d
!     MODULE PROCEDURE SpAMM_Add_Scalar_2_Tree_2d_full
!     MODULE PROCEDURE SpAMM_Add_Scalar_2_Tree_2d
     !--------------------------------------------------------
!     MODULE PROCEDURE SpAMM_Add_QuTree_2_QuTree_InPlace
!     MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_InPlace
!     MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_RePlace
!     MODULE PROCEDURE SpAMM_Add_Identity_2_QuTree_InPlace
     !--------------------------------------------------------
!     module procedure spamm_add_identity_to_matrix_2nd_order
!     module procedure spamm_add_2nd_order_to_2nd_order
  END INTERFACE Add
CONTAINS

  FUNCTION SpAMM_Add_Tree2d_2_Tree2d (A, B, alpha, beta, C) RESULT(D)   
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT)           :: A, B
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT), OPTIONAL :: C
    TYPE(SpAMM_Tree_2d), POINTER                          :: D
    real(spamm_kind), intent(in), optional                :: alpha, beta
    REAL(SpAMM_KIND)                                      :: Local_Alpha,Local_Beta
    !
    D=>NULL()
    !
    ! hmmm... this needs some thought ...
    !if(.not. associated(A))RETURN
    !if(.not. associated(B))RETURN
    !
    IF(PRESENT(Alpha))THEN
       Local_Alpha=Alpha
    ELSE
       Local_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
       Local_Beta=Beta
    ELSE
       Local_Beta=SpAMM_One
    ENDIF

    IF(PRESENT(C))THEN
       !  C=>Add(A,B,a,C): in place on C, C can be A, B or other.
       IF(ASSOCIATED(B,C))THEN
          ! B=>a*A+b*B
          CALL SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur(B, A, Local_beta, Local_alpha)
          D=>B
       ELSEIF(ASSOCIATED(A,C))THEN
          ! A=>a*A+b*B
          CALL SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur(A, B, Local_alpha, Local_beta)
          D=>A
       ELSE
          ! C=>a*A+B
          CALL SpAMM_Add_Tree2d_2_Tree2d_Recur(C, A, B, Local_alpha, Local_beta)
          D=>C
       ENDIF
    ELSE
       ! D=>a*A+B
       !       How to new here also needs some thought ..
       !       CALL New(D, A%decoration )
       CALL SpAMM_Add_Tree2d_2_Tree2d_Recur(D, A, B, Local_alpha, Local_beta)       
    ENDIF
  END FUNCTION SpAMM_Add_Tree2d_2_Tree2d

  RECURSIVE SUBROUTINE SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur(A, B, alpha, beta)
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT) :: A
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(IN)    :: B
    REAL(SpAMM_KIND),             INTENT(IN)    :: alpha, beta
    logical                                     :: TA, TB
    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)
    IF(.NOT.TA.AND.TB)THEN
       CALL SpAMM_Copy_Tree2d_2_Tree2d(B, A)
       CALL SpAMM_Multiply_Tree2d_x_Scalar_recur(A, beta)
    ELSEIF(TA .AND. .NOT.TB) THEN
       ! Multiply A tree with alpha.
       CALL SpAMM_Multiply_Tree2d_x_Scalar_recur(A, alpha)
    ELSEIF(TA.AND.TB)THEN
       IF(ALLOCATED(A%Blok).AND.ALLOCATED(B%Blok))THEN
          A%Blok=alpha*A%Blok+beta*B%Blok
       ELSE
          CALL SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur(A%child_00, B%child_00, alpha, beta)
          CALL SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur(A%child_01, B%child_01, alpha, beta)
          CALL SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur(A%child_11, B%child_11, alpha, beta)
       ENDIF
       CALL Decorate(A)
    ENDIF
  END SUBROUTINE SpAMM_Add_Tree2d_2_Tree2d_InPlace_Recur

  RECURSIVE SUBROUTINE SpAMM_Add_Tree2d_2_Tree2d_Recur(A, B, C, alpha, beta)
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(IN)    :: A,B
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT) :: C
    REAL(SpAMM_KIND)                            :: alpha, beta
    logical                                     :: TA, TB
    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)
    IF(.NOT.TA.AND.TB)THEN
       CALL SpAMM_Copy_Tree2d_2_Tree2d(B, A)
       CALL SpAMM_Multiply_Tree2d_x_Scalar_recur(A, beta)
    ELSEIF(TA .AND. .NOT.TB) THEN
       ! Multiply A tree with alpha.
       CALL SpAMM_Multiply_Tree2d_x_Scalar_recur(A, alpha)
    ELSEIF(TA.AND.TB)THEN
       IF(ALLOCATED(A%Blok).AND.ALLOCATED(B%Blok))THEN
          A%Blok=alpha*A%Blok+beta*B%Blok
       ELSE
          CALL SpAMM_Add_Tree2d_2_Tree2d_Recur(A%child_00, B%child_00, C%child_00, alpha, beta)
          CALL SpAMM_Add_Tree2d_2_Tree2d_Recur(A%child_01, B%child_01, C%child_01, alpha, beta)
          CALL SpAMM_Add_Tree2d_2_Tree2d_Recur(A%child_11, B%child_11, C%child_11, alpha, beta)
       ENDIF
       CALL Decorate(C)
    ENDIF
  END SUBROUTINE SpAMM_Add_Tree2d_2_Tree2d_Recur

  FUNCTION SpAMM_Add_Tree1d_2_Tree1d (A, B, alpha, beta, C) RESULT(D)   
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT)           :: A, B
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT), OPTIONAL :: C
    TYPE(SpAMM_Tree_1d), POINTER                          :: D
    real(spamm_kind), intent(in), optional                :: alpha, beta
    REAL(SpAMM_KIND)                                      :: Local_Alpha,Local_Beta
    !
    D=>NULL()
    !
    ! hmmm... this needs some thought ...
    !if(.not. associated(A))RETURN
    !if(.not. associated(B))RETURN
    !
    IF(PRESENT(Alpha))THEN
       Local_Alpha=Alpha
    ELSE
       Local_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
       Local_Beta=Beta
    ELSE
       Local_Beta=SpAMM_One
    ENDIF

    IF(PRESENT(C))THEN
       !  C=>Add(A,B,a,C): in place on C, C can be A, B or other.
       IF(ASSOCIATED(B,C))THEN
          ! B=>a*A+b*B
          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(B, A, Local_beta, Local_alpha)
          D=>B
       ELSEIF(ASSOCIATED(A,C))THEN
          ! A=>a*A+b*B
          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A, B, Local_alpha, Local_beta)
          D=>A
       ELSE
          ! C=>a*A+B
          CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(C, A, B, Local_alpha, Local_beta)
          D=>C
       ENDIF
    ELSE
       ! D=>a*A+B
       !       CALL New(D, A%decoration )
       CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(D, A, B, Local_alpha, Local_beta)       
    ENDIF
  END FUNCTION SpAMM_Add_Tree1d_2_Tree1d

  RECURSIVE SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A, B, alpha, beta)
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT) :: A
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(IN)    :: B
    REAL(SpAMM_KIND),             INTENT(IN)    :: alpha, beta
    logical                                     :: TA, TB
    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)
    IF(.NOT.TA.AND.TB)THEN
       CALL SpAMM_Copy_Tree1d_2_Tree1d(B, A)
       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, beta)
    ELSEIF(TA .AND. .NOT.TB) THEN
       ! Multiply A tree with alpha.
       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, alpha)
    ELSEIF(TA.AND.TB)THEN
       IF(ALLOCATED(A%Blok).AND.ALLOCATED(B%Blok))THEN
          A%Blok=alpha*A%Blok+beta*B%Blok
       ELSE
          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A%child_0, B%child_0, alpha, beta)
          CALL SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur(A%child_1, B%child_1, alpha, beta)
       ENDIF
       CALL ReDecorate(A)
    ENDIF
  END SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_InPlace_Recur

  RECURSIVE SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_Recur(A, B, C, alpha, beta)
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(IN)    :: A,B
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT) :: C
    REAL(SpAMM_KIND)                            :: alpha, beta
    logical                                     :: TA, TB
    TA=ASSOCIATED(A)
    TB=ASSOCIATED(B)
    IF(.NOT.TA.AND.TB)THEN
       CALL SpAMM_Copy_Tree1d_2_Tree1d(B, A)
       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, beta)
    ELSEIF(TA .AND. .NOT.TB) THEN
       ! Multiply A tree with alpha.
       CALL SpAMM_Multiply_Tree1d_x_Scalar_recur(A, alpha)
    ELSEIF(TA.AND.TB)THEN
       IF(ALLOCATED(A%Blok).AND.ALLOCATED(B%Blok))THEN
          A%Blok=alpha*A%Blok+beta*B%Blok
       ELSE
          CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(A%child_0, B%child_0, C%child_0, alpha, beta)
          CALL SpAMM_Add_Tree1d_2_Tree1d_Recur(A%child_1, B%child_1, C%child_1, alpha, beta)
       ENDIF
       CALL ReDecorate(C)
    ENDIF
  END SUBROUTINE SpAMM_Add_Tree1d_2_Tree1d_Recur


!  INCLUDE 'spamm_algebra_add_old.F90'

END module spamm_algebra_add
