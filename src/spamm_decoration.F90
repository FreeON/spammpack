! The "long weekend hack".  Going for a functional implementation
! that fits tightly in the edit window (mc,2/2015).
!
module spamm_decoration

  use  spamm_structures
  !
  implicit none

!  PRIVATE
  ! some vanilla ... 
!  PUBLIC :: Decorate
!  INTERFACE Decorate
!     MODULE PROCEDURE SpAMM_redecorate_1d
!     MODULE PROCEDURE SpAMM_redecorate_2d
!  END INTERFACE Decorate

CONTAINS

  ! Protocols for the uppward, child -> parent merge of SpAMM Tree Decorations, STDECs: 
  ! assuming bounding box has been set first in a dowards pass, we should now just be in 
  ! the backwards accumulation phase for upwards merge with these data structures

  ! - - - - - - - - - - - - - - - - -1d, 1d, 1d - - - - - - - - - - - - - - - - - - 

  ! uppwards redecoration ...
  SUBROUTINE SpAMM_redecorate_tree_1d(a)

    TYPE(SpAMM_Tree_1d), POINTER, INTENT(INOUT) :: a    


    if(a%frill%leaf)then  ! at a leaf?
       a%frill%Norm2=SUM(a%chunk**2) 
       a%frill%Non0s=SBS
       ! application has to fill in the %flops at this level
       RETURN
    ELSE ! init this level
       a%frill%Norm2=SpAMM_Zero
       a%frill%Non0s=SpAMM_Zero
       a%frill%FlOps=SpAMM_Zero
    ENDIF

!    WRITE(*,*)'-----------------------------------------------'
!    WRITE(*,*)' bb = ',a%frill%bndbx
!    WRITE(*,*)'walking back up ...', a%frill%norm2

    ! walk back up one level with each decoration ...
    IF(ASSOCIATED( a%child_0 )) CALL SpAMM_merge_1d( a%frill, a%child_0%frill )
    IF(ASSOCIATED( a%child_1 )) CALL SpAMM_merge_1d( a%frill, a%child_1%frill )

!    if(a%child_0%frill%norm2<0d0)stop '0'
!    if(a%child_1%frill%norm2<0d0)then
!       write(*,*)' asociated ',ASSOCIATED( a%child_1 ),a%child_1%frill%leaf
!       write(*,*)a%child_1%frill%bndbx
!       write(*,*)a%child_1%frill%norm2
!       stop '1'
!    endif

!    WRITE(*,*)'done now ', a%frill%norm2


  END SUBROUTINE SpAMM_redecorate_tree_1d
 
  ! 1d decoration merge:
  SUBROUTINE SpAMM_merge_1d(a,b)

    TYPE(SpAMM_decoration_1d), INTENT(INOUT) :: a
    TYPE(SpAMM_decoration_1d), INTENT(IN)    :: b

    a%Norm2=a%Norm2+b%Norm2
        
!    write(*,*)' parent bb = ',a%bndbx
!    write(*,*)' child  bb = ',b%bndbx
!    write(*,*)' child norm= ',b%norm2

    a%Non0s=a%Non0s+b%Non0s
    a%FlOps=a%FlOps+b%FlOps

  END SUBROUTINE SpAMM_merge_1d

  !  d_2d |cpy> a_2d (can be used top down or bottom up...)
  SUBROUTINE SpAMM_decoration_1d_copy_decoration_1d(d,a)

    type(SpAMM_decoration_1d) :: d
    type(SpAMM_decoration_1d) :: a

    d%leaf =a%leaf
    d%norm2=a%norm2
    d%non0s=a%non0s
    d%flops=a%flops
    d%ndimn=a%ndimn
    d%bndbx=a%bndbx

  END SUBROUTINE SpAMM_decoration_1d_copy_decoration_1d


  ! - - - - - - - - - - - - - - - - -2d, 2d, 2d - - - - - - - - - - - - - - - - - - 

  !  d_2d |cpy> a_2d (can be used top down or bottom up...)
  SUBROUTINE SpAMM_decoration_2d_copy_decoration_2d(d,a)

    type(SpAMM_decoration_2d) :: d
    type(SpAMM_decoration_2d) :: a

    d%leaf =a%leaf
    d%norm2=a%norm2
    d%non0s=a%non0s
    d%flops=a%flops
    d%ndimn=a%ndimn
    d%bndbx=a%bndbx

  END SUBROUTINE SpAMM_decoration_2d_copy_decoration_2d

  ! uppwards redecoration ...
  SUBROUTINE SpAMM_redecorate_tree_2d_symm(a)

    TYPE(SpAMM_tree_2d_symm), POINTER,  INTENT(INOUT) :: a    

    if(a%frill%leaf)then  ! at a leaf?
       a%frill%Norm2=SUM(a%chunk**2) 
       a%frill%Non0s=SBS2
       ! application has to fill in the %flops at this level
       RETURN
    ELSE ! init this level
       a%frill%Norm2=SpAMM_Zero
       a%frill%Non0s=SpAMM_Zero
       a%frill%FlOps=SpAMM_Zero
    ENDIF

    ! walk back up one level with each decoration ...
    IF(ASSOCIATED(A%child_00))CALL SpAMM_merge_decoration_2d(A%frill, A%child_00%frill)
    IF(ASSOCIATED(A%child_01))CALL SpAMM_merge_decoration_2d(A%frill, A%child_01%frill)
    IF(ASSOCIATED(A%child_11))CALL SpAMM_merge_decoration_2d(A%frill, A%child_11%frill)

  END SUBROUTINE SpAMM_redecorate_tree_2d_symm

  ! the 2d decoration merge:
  SUBROUTINE SpAMM_merge_decoration_2d(a,b)

    TYPE(SpAMM_decoration_2d), INTENT(INOUT) :: a
    TYPE(SpAMM_decoration_2d), INTENT(IN)    :: b

    a%Norm2=a%Norm2+b%Norm2
    a%Non0s=a%Non0s+b%Non0s
    a%FlOps=a%FlOps+b%FlOps

  END SUBROUTINE SpAMM_merge_decoration_2d

END module spamm_decoration

