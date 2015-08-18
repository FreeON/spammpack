module spamm_decorate
  use spamm_globals
  use spamm_types
  implicit none
  PRIVATE
  PUBLIC :: Decorate
  !> Interface for additions operations between different SpAMM types.
  INTERFACE Decorate
     MODULE PROCEDURE SpAMM_Decorate_Tree1d
     MODULE PROCEDURE SpAMM_Decorate_Tree2d
  END INTERFACE Decorate
CONTAINS

  SUBROUTINE SpAMM_Decorate_Tree2d(A)
    TYPE(SpAMM_Tree_2d), POINTER, INTENT(INOUT)   :: A    
    IF(ALLOCATED(A%Blok))THEN
       A%decoration%norm2=SUM(A%Blok**2)
       A%decoration%number_nonzeros=SIZE(A%Blok,1)*SIZE(A%Blok,2)
       RETURN
    ENDIF
    IF(ASSOCIATED(A%child_00))CALL DecorateUp_2d(A%decoration,A%child_00%decoration)
    IF(ASSOCIATED(A%child_01))CALL DecorateUp_2d(A%decoration,A%child_01%decoration)
    IF(ASSOCIATED(A%child_11))CALL DecorateUp_2d(A%decoration,A%child_11%decoration)
  END SUBROUTINE SpAMM_Decorate_Tree2d

  SUBROUTINE SpAMM_Decorate_Tree1d(A)
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT)   :: A    
    IF(ALLOCATED(A%Vect))THEN
       A%decoration%norm2=SUM(A%Vect**2)
       A%decoration%number_nonzeros=SIZE(A%Vect,1)
       RETURN
    ENDIF
    IF(ASSOCIATED(A%child_0))CALL DecorateUp_1d(A%decoration,A%child_0%decoration)
    IF(ASSOCIATED(A%child_1))CALL DecorateUp_1d(A%decoration,A%child_1%decoration)
  END SUBROUTINE SpAMM_Decorate_Tree1d

  SUBROUTINE DecorateUp_2d(A,B)
    TYPE(decoration_2d), INTENT(INOUT) :: A
    TYPE(decoration_2d), INTENT(IN)    :: B
    A%bb(1,0)=MIN(A%bb(1,0),B%bb(1,0))
    A%bb(2,0)=MAX(A%bb(2,0),B%bb(2,0))
    A%bb(1,1)=MIN(A%bb(1,1),B%bb(1,1))
    A%bb(2,1)=MAX(A%bb(2,1),B%bb(2,1))
    A%norm2=A%norm2+B%norm2
    A%number_nonzeros=A%number_nonzeros+B%number_nonzeros    
  END SUBROUTINE DecorateUp_2d

  SUBROUTINE DecorateUp_1d(A,B)
    TYPE(decoration_1d), INTENT(INOUT) :: A
    TYPE(decoration_1d), INTENT(IN)    :: B
    A%bb(1)=MIN(A%bb(1),B%bb(1))
    A%bb(2)=MAX(A%bb(2),B%bb(2))
    A%norm2=A%norm2+B%norm2
    A%number_nonzeros=A%number_nonzeros+B%number_nonzeros    
  END SUBROUTINE DecorateUp_1d

END module spamm_decorate
