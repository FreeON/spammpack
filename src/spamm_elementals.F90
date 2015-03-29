module spamm_elementals

   use spamm_structures
   use spamm_xstructors
   use spamm_decoration

  implicit none

  ! This code contains protocols for solo/unit opperations on structures, copy, scalar multiply etc. 
CONTAINS


  recursive subroutine SpAMM_set_identity_2d_symm_recur (A)

    type(SpAMM_tree_2d_symm), pointer  :: A
    integer :: i

    if(.not. associated(A))then
       return
    else

       if(a%frill%leaf)then
          
          a%frill%init=.FALSE. 
          a%chunk=SpAMM_zero
          DO I=1,a%frill%bndbx(1,1)-a%frill%bndbx(0,1)+1
             a%chunk(I,I)=SpAMM_one
          ENDDO

       else

          CALL SpAMM_set_identity_2d_symm_recur(SpAMM_construct_tree_2d_symm_00(a))
          call SpAMM_destruct_tree_2d_symm_recur (a%child_01)
          call SpAMM_destruct_tree_2d_symm_recur (a%child_10)
          CALL SpAMM_set_identity_2d_symm_recur(SpAMM_construct_tree_2d_symm_11(a))

       endif

    endif

    CALL SpAMM_redecorate_tree_2d_symm(a)

  end subroutine SpAMM_set_identity_2d_symm_recur


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
       CALL SpAMM_print_tree_2d_symm_recur (a%child_10) 
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
          absmax = max(absmax, SpAMM_absmax_tree_2d_symm_recur(a%child_10) )
          absmax = max(absmax, SpAMM_absmax_tree_2d_symm_recur(a%child_11) )

       endif

    endif

  end function SpAMM_absmax_tree_2d_symm_recur

  recursive function SpAMM_twist_2d_symm_recur (A, at_o) result(twist)

    type(SpAMM_tree_2d_symm), pointer,           intent(in)   :: A
    type(SpAMM_tree_2d_symm), pointer, optional, intent(in)   :: At_o
    real(SpAMM_KIND)                                :: twist, twist_00, twist_11, twist_01, twist_10

    twist  = 0

    if(.not. associated(a))return

       if(a%frill%leaf)then

          IF(PRESENT(at_o))THEN
             
             twist=SUM( (a%chunk(1:SBS,1:SBS)-TRANSPOSE(at_o%chunk(1:SBS,1:SBS)))**2 )

          ELSE

             twist=SUM( (a%chunk(1:SBS,1:SBS)-TRANSPOSE(   a%chunk(1:SBS,1:SBS)))**2 )

          ENDIF

       else

          twist_00=SpAMM_twist_2d_symm_recur( a%child_00 ) 
          twist_11=SpAMM_twist_2d_symm_recur( a%child_11 )

          twist_01=SpAMM_twist_2d_symm_recur( a%child_01, a%child_10 )
          twist_10=SpAMM_twist_2d_symm_recur( a%child_10, a%child_01 )

          twist = twist_00 + twist_11 + twist_01 + twist_10

       endif

  end function SpAMM_twist_2d_symm_recur

END module spamm_elementals
