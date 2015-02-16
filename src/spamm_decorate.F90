! The "long weekend hack".  Going for a functional implementation
! that fits tightly in the edit window (mc,2/2015).
!
! ********* The SpAMM Garnish (SGARN) tree enrichments ...
module spamm_decorate

  use  spamm_globals
  use  spamm_types
  !
  implicit none
  PRIVATE

  ! some vanilla ... 
  PUBLIC :: Decorate
  INTERFACE Decorate
     MODULE PROCEDURE SpAMM_redecorate_1d
     MODULE PROCEDURE SpAMM_redecorate_2d
  END INTERFACE Decorate

CONTAINS

  ! Protocols for the uppward, child->parent merge of the SGARN ... 

  ! - - - - - - - - - - - - - - - - -1d, 1d, 1d - - - - - - - - - - - - - - - - - - 
  ! uppwards redecoration ...
  SUBROUTINE SpAMM_redecorate_1d(a)
    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT) :: a    
    ! at a leaf ...
    IF(ALLOCATED(a%chunk))THEN 
       a%norm2=SUM(a%chunk**2) 
       a%non0s=SIZE(a%chunk)
       RETURN 
    ENDIF
    ! re-init this level
    a%frill%Norm2=0; a%frill%Non0s=0; a%frill%BndBx=0
    ! walk back up one level with each decoration ...
    IF(ASSOCIATED( a%child_0 )) CALL SpAMM_merge_1d( a%frill, a%child_0%frill )
    IF(ASSOCIATED( a%child_1 )) CALL SpAMM_merge_1d( a%frill, a%child_1%frill )
  END SUBROUTINE SpAMM_redecorate_1d   
 
  ! 1d decoration merge:
  SUBROUTINE SpAMM_merge_1d(a,b)
    TYPE(decoration_1d), INTENT(INOUT) :: a
    TYPE(decoration_1d), INTENT(IN)    :: b
    a%Norm2=a%Norm2+a%Norm2
    a%Non0s=a%Non0s+a%Non0s
    !    a%BndBx(0)=MIN(a%BndBx(0),b%BndBx(0))
    !    a%BndBx(1)=MAX(a%BndBx(1),b%BndBx(1))
  END SUBROUTINE SpAMM_merge_1d
  !
  ! - - - - - - - - - - - - - - - - -2d, 2d, 2d - - - - - - - - - - - - - - - - - - 
  ! uppwards redecoration ...
  SUBROUTINE SpAMM_redecorate_2d_symm(a)

    TYPE(SpAMM_Tree_2d_symm), POINTER, INTENT(INOUT)   :: a    

    ! at a leaf?
    IF(ALLOCATED(a%chunk))THEN 
       a%norm2=SUM(a%chunk*2) 
       a%Non0s=SIZE(a%chunk,1)*SIZE(a%chunk,2)
       RETURN
    ENDIF

    ! re-init at this level
    a%frill%Norm2=SpAMM_Zero
    a%frill%Non0s=SpAMM_Zero
    a%frill%FlOps=SpAMM_Zero
    !    a%frill%BndBx=SpAMM_Zero

    ! walk back up one level with each decoration ...
    IF(ASSOCIATED(A%child_00))CALL SpAMM_merge_decoration_2d(A%frill, A%child_00%frill)
    IF(ASSOCIATED(A%child_01))CALL SpAMM_merge_decoration_2d(A%frill, A%child_01%frill)
    IF(ASSOCIATED(A%child_11))CALL SpAMM_merge_decoration_2d(A%frill, A%child_11%frill)

  END SUBROUTINE SpAMM_redecorate_2d_symm

  ! the 2d decoration merge:
  SUBROUTINE SpAMM_merge_decoration_2d(a,b)

    TYPE(decoration_2d), INTENT(INOUT) :: a
    TYPE(decoration_2d), INTENT(IN)    :: b
    ! assuming bounding box has been set first in a dowards pass ...
    ! we should now just be in the backwards accumulation phase at this point ...
    a%Norm2=a%Norm2+b%Norm2
    a%Non0s=a%Non0s+b%Non0s
    a%FlOps=a%FlOps+b%FlOps
  END SUBROUTINE SpAMM_merge_decoration_2d

    !    a%BndBx(0,1)=MIN(a%BndBx(0,1),b%BndBx(0,1))
    !    a%BndBx(0,2)=MAX(a%BndBx(0,2),b%BndBx(0,2))
    !    a%BndBx(1,1)=MIN(a%BndBx(1,1),b%BndBx(1,1))
    !    a%BndBx(1,2)=MAX(a%BndBx(1,2),b%BndBx(1,2))

END module spamm_decorate
