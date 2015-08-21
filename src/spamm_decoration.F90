! The "long weekend hack".  Going for a functional implementation
! that fits tightly in the edit window (mc,2/2015).
!
module spamm_decoration

  use  spamm_structures

  !  use  spamm_xstructors

  !
  implicit none

contains

  ! Protocols for the uppward, child -> parent merge of SpAMM Tree Decorations, STDECs:
  ! assuming bounding box has been set first in a dowards pass, we should now just be in
  ! the backwards accumulation phase for upwards merge with these data structures

  ! - - - - - - - - - - - - - - - - -1d, 1d, 1d - - - - - - - - - - - - - - - - - -

  ! uppwards redecoration ...
  subroutine SpAMM_redecorate_tree_1d(a)

    type(SpAMM_Tree_1d), pointer :: a

    if(.not.associated(a))return

    if(a%frill%leaf)then  ! at a leaf?
       a%frill%Norm2=sum(a%chunk**2)
       a%frill%Non0s=SBS
       ! application has to fill in the %flops at this level
       return
    else ! init this level
       a%frill%Norm2=SpAMM_Zero
       a%frill%Non0s=SpAMM_Zero
       a%frill%FlOps=SpAMM_Zero
    end if

    ! walk back up one level with each decoration ...
    call SpAMM_merge_1d( a, a%child_0 )
    call SpAMM_merge_1d( a, a%child_1 )

  end subroutine SpAMM_redecorate_tree_1d

  ! 1d decoration merge:
  subroutine SpAMM_merge_1d(a,b)

    type(SpAMM_Tree_1d), pointer :: a,b

    ! nothing to be seen here ...
    if(.not.associated(b)) return

    ! this is unused, passed in data, don't accumulate
    if(b%frill%init ) return

    a%frill%Init =a%frill%Init .and. b%frill%Init
    a%frill%Norm2=a%frill%Norm2+b%frill%Norm2
    a%frill%Non0s=a%frill%Non0s+b%frill%Non0s
    a%frill%FlOps=a%frill%FlOps+b%frill%FlOps

  end subroutine SpAMM_merge_1d

  !  d_2d |cpy> a_2d (can be used top down or bottom up...)

  ! uppwards redecoration ...
  subroutine SpAMM_redecorate_tree_2d_symm(a)

    type(SpAMM_tree_2d_symm), pointer :: a

    if(.not.associated(a)) return

    if(a%frill%leaf)then  ! at a leaf?
       a%frill%Norm2=sum(a%chunk**2)
       a%frill%Non0s=SBS2
       ! application has to fill in the %flops at this level
       return
    else ! init this level
       a%frill%Norm2=SpAMM_Zero
       a%frill%Non0s=SpAMM_Zero
       a%frill%FlOpS=SpAMM_Zero
    end if

    call SpAMM_merge_decoration_2d(a,a%child_00)
    !    if(associated(a%child_00)) &
    !    WRITE(*,*)' a norm00 = ',a%frill%norm2,a%child_00%frill%norm2

    call SpAMM_merge_decoration_2d(a,a%child_11)
    !    if(associated(a%child_11)) &
    !    WRITE(*,*)' a norm11 = ',a%frill%norm2,a%child_11%frill%norm2

    call SpAMM_merge_decoration_2d(a,a%child_01)
    !    if(associated(a%child_01)) &
    !    WRITE(*,*)' a norm01 = ',a%frill%norm2,a%child_01%frill%norm2

    call SpAMM_merge_decoration_2d(a,a%child_10)
    !    if(associated(a%child_10)) &
    !    WRITE(*,*)' a norm10 = ',a%frill%norm2,a%child_10%frill%norm2


  end subroutine SpAMM_redecorate_tree_2d_symm

  ! the 2d decoration merge:
  subroutine SpAMM_merge_decoration_2d(a,b)

    type(SpAMM_tree_2d_symm), pointer :: a,b

    ! nothing to be seen here ...
    if(.not.associated(b)) return

    ! this is unused, passed in data, don't accumulate
    if(b%frill%init ) return

    !    write(*,*)' flops = ',a%frill%FlOps,b%frill%FlOps

    a%frill%Init =a%frill%Init .and. b%frill%Init
    a%frill%Norm2=a%frill%Norm2+b%frill%Norm2
    a%frill%Non0s=a%frill%Non0s+b%frill%Non0s
    a%frill%FlOps=a%frill%FlOps+b%frill%FlOps


  end subroutine SpAMM_merge_decoration_2d

end module spamm_decoration
