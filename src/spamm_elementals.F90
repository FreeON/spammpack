module spamm_elementals

   use spamm_structures
   use spamm_xstructors
   use spamm_decoration

  implicit none

  ! This code contains protocols for solo opperations on structures, copy, scalar multiply etc. 
CONTAINS

  RECURSIVE FUNCTION SpAMM_trace_tree_2d_symm_recur(a) RESULT(tr)
    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN) :: a
    REAL(SpAMM_KIND)                              :: tr, tr00,tr11
    integer                                       :: i

    tr=SpAMM_Zero
    IF(.NOT.ASSOCIATED(a))RETURN

    IF(a%frill%leaf)THEN ! Leaf condition ? 

       do i = 1, a%frill%bndbx(1,1) - a%frill%bndbx(0,1)+1
         tr = tr + a%chunk(i,i)
       enddo

    ELSE

       tr00=SpAMM_trace_tree_2d_symm_recur(a%child_00) 
       tr11=SpAMM_trace_tree_2d_symm_recur(a%child_11)
       tr=tr00+tr11

    ENDIF
  END FUNCTION SpAMM_trace_tree_2d_symm_recur

  recursive subroutine SpAMM_print_tree_2d_symm_recur (A) 

    type(SpAMM_tree_2d_symm), pointer, intent(in)   :: A
    integer :: i,k

    if(.not.associated(A))return

    if(a%frill%leaf.and.a%frill%norm2>1d-3)then

       WRITE(*,33)a%frill%bndbx(0:1,1),a%frill%bndbx(0:1,2)
33     FORMAT(' [',I4,", ",I4,"]x[",I4,", ",I4,"]")
       DO I=1,SBS
          WRITE(*,44)(a%chunk(k,I),k=1,SBS)
       ENDDO
44     FORMAT(20(E10.4,", "))

    else

       CALL SpAMM_print_tree_2d_symm_recur (a%child_00) 
       CALL SpAMM_print_tree_2d_symm_recur (a%child_01) 
       CALL SpAMM_print_tree_2d_symm_recur (a%child_11) 

    endif

  end subroutine SpAMM_print_tree_2d_symm_recur



  recursive subroutine SpAMM_print_tree_1d_recur (A) 

    type(SpAMM_tree_1d), pointer, intent(in)   :: A

    if(.not.associated(A))return
    
    if(a%frill%leaf)then
       if(sum(a%chunk**2)<1d-6)return
       WRITE(*,33)a%frill%bndbx(0:1), a%chunk
33     FORMAT(' [',I4,"-",I4,"]: ",20(E10.4,", ") )

    else

       CALL SpAMM_print_tree_1d_recur (a%child_0) 
       CALL SpAMM_print_tree_1d_recur (a%child_1) 

    endif

  end subroutine SpAMM_print_tree_1d_recur

  recursive function SpAMM_absmax_tree_2d_symm_recur (A) result(absmax)

    type(SpAMM_tree_2d_symm), pointer, intent(in)   :: A
    real(SpAMM_KIND)                                :: absmax
    integer, dimension(:,:),  pointer               :: bb

    absmax = 0
    if(.not. associated(A))then
       return
    else

       bb => A%frill%bndbx 
       if(bb(1,1)-bb(0,1)==SBS)then

          absmax = maxval(abs(a%chunk))

       else

          absmax = max(absmax, SpAMM_absmax_tree_2d_symm_recur(a%child_00) )
          absmax = max(absmax, SpAMM_absmax_tree_2d_symm_recur(a%child_01) )
          absmax = max(absmax, SpAMM_absmax_tree_2d_symm_recur(a%child_11) )

       endif

    endif

  end function SpAMM_absmax_tree_2d_symm_recur




END module spamm_elementals
