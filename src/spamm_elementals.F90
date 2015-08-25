module spamm_elementals

   use spamm_structures
   use spamm_xstructors
   use spamm_decoration

  implicit none

  ! This code contains protocols for solo/unit opperations on structures, copy, scalar multiply etc.
CONTAINS


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (wrapper)
  FUNCTION SpAMM_set_identity_2d_symm (MN, in_O ) RESULT(d)

    INTEGER, DIMENSION(2), INTENT(IN) :: MN
    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL :: In_O
    TYPE(SpAMM_tree_2d_symm), POINTER           :: D
    INTEGER                                     :: Depth

    ! figure the starting conditions ...
    if(present(in_O))then
       d => in_O
    else
       d => NULL()
    endif

    if(.not.associated(d))then
       ! instantiate a tree if no passed allocation
       d => SpAMM_new_top_tree_2d_symm(MN)
    endif

    ! set passed data for initialization
    CALL SpAMM_flip(d)

    depth=0
    CALL SpAMM_set_identity_2d_symm_recur (d, depth)

    ! prune unused nodes ...
    CALL SpAMM_prune(d)

  END FUNCTION SpAMM_set_identity_2d_symm




  recursive subroutine SpAMM_set_identity_2d_symm_recur (A, depth)

    type(SpAMM_tree_2d_symm), pointer  :: A
    integer :: depth, i

    if(.not. associated(A))then
       return
    else

       if(a%frill%leaf)then

          a%frill%needs_initialization=.FALSE.
          a%chunk=0
          DO I=1,a%frill%bndbx(1,1)-a%frill%bndbx(0,1)+1
             a%chunk(I,I)=SpAMM_one
          ENDDO

       else

          CALL SpAMM_set_identity_2d_symm_recur(SpAMM_construct_tree_2d_symm_00(a), depth+1)
          CALL SpAMM_set_identity_2d_symm_recur(SpAMM_construct_tree_2d_symm_11(a), depth+1)

       endif

    endif

    CALL SpAMM_redecorate_tree_2d_symm(a)

  end subroutine SpAMM_set_identity_2d_symm_recur


  RECURSIVE FUNCTION SpAMM_trace_tree_2d_symm_recur(a) RESULT(tr)
    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN) :: a
    REAL(SPAMM_KIND)                              :: tr, tr00,tr11
    integer                                       :: i

    tr=0
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

       WRITE(*,33)a%frill%bndbx(0:1,1),a%frill%bndbx(0:1,2),a%frill%norm2
33     FORMAT(' [',I4,", ",I4,"]x[",I4,", ",I4,"] = ",e12.6)
!       DO I=1,SBS
!          WRITE(*,44)(a%chunk(k,I),k=1,SBS)
!       ENDDO
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
    real(SPAMM_KIND)                                :: absmax
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


  recursive function SpAMM_twist_tree_2d_symm_double_recur (A, AT) result(twist)

    type(SpAMM_tree_2d_symm), pointer  :: A,AT
    real(SPAMM_KIND)                                :: twist, twist_00, twist_11, twist_01, twist_10

    twist  = 0

    if(.not. associated(a))return
    if(.not. associated(at))return

    if(a%frill%leaf)then

       twist=SUM( (a%chunk(1:SBS,1:SBS)-TRANSPOSE(at%chunk(1:SBS,1:SBS)))**2 )

    ELSE

       twist_01=SpAMM_twist_tree_2d_symm_double_recur ( a%child_01, a%child_10 )
       twist_10=SpAMM_twist_tree_2d_symm_double_recur ( a%child_10, a%child_01 )

!       twist_00=SpAMM_twist_2d_symm_single_recur( a%child_00 )
!       twist_11=SpAMM_twist_2d_symm_single_recur( a%child_11 )

       twist = twist_01 + twist_10 !twist_00 + twist_11 +

    ENDIF

  end function SpAMM_twist_tree_2d_symm_double_recur

  function SpAMM_twist_tree_2d_symm(A) result(twist)

    type(SpAMM_tree_2d_symm), pointer  :: A
    real(SPAMM_KIND)                                :: twist, twist_00, twist_11, twist_01, twist_10

    twist  = 0

    if(.not. associated(a))return

       twist_01=SpAMM_twist_tree_2d_symm_double_recur ( a%child_01, a%child_10 )
       twist_10=SpAMM_twist_tree_2d_symm_double_recur ( a%child_10, a%child_01 )

!       twist_00=SpAMM_twist_2d_symm_single_recur( a%child_00 )
!       twist_11=SpAMM_twist_2d_symm_single_recur( a%child_11 )

       twist=twist_01+twist_10

       twist=SQRT(twist)

     end function SpAMM_twist_tree_2d_symm

  recursive subroutine spamm_tree_print_leaves_2d_symm_recur(A, file_unit)

    type(spamm_tree_2d_symm), pointer, intent(in)   :: A
    integer, intent(in) :: file_unit

    integer :: i,k

    if(.not.associated(A)) return

    if(A%frill%leaf) then
       write(file_unit, "(2ES20.10,2I8,ES20.10)") &
            A%frill%bndbx(0,1)+real(A%frill%bndbx(1,1)-A%frill%bndbx(0,1)+1)/2., &
            A%frill%bndbx(0,2)+real(A%frill%bndbx(1,2)-A%frill%bndbx(0,2)+1)/2., &
            A%frill%bndbx(1,1)-A%frill%bndbx(0,1)+1, &
            A%frill%bndbx(1,2)-A%frill%bndbx(0,2)+1, &
            LOG10(SQRT(A%frill%norm2))
    else
       call spamm_tree_print_leaves_2d_symm_recur(A%child_00, file_unit)
       call spamm_tree_print_leaves_2d_symm_recur(A%child_01, file_unit)
       call spamm_tree_print_leaves_2d_symm_recur(A%child_10, file_unit)
       call spamm_tree_print_leaves_2d_symm_recur(A%child_11, file_unit)
    endif

  end subroutine spamm_tree_print_leaves_2d_symm_recur

  !> Print bounding boxes and norms of the leave nodes. Use for
  !> visualization.
  !!
  !! \param A The matrix.
  !! \param file_unit The optional unit of the output file.
  subroutine spamm_tree_print_leaves_2d_symm(A, file_unit)

    type(spamm_tree_2d_symm), pointer, intent(in) :: A
    integer, intent(in) :: file_unit

    write(file_unit, "(I8)") SPAMM_CHUNK_SIZE

    call spamm_tree_print_leaves_2d_symm_recur(A%child_00, file_unit)
    call spamm_tree_print_leaves_2d_symm_recur(A%child_01, file_unit)
    call spamm_tree_print_leaves_2d_symm_recur(A%child_10, file_unit)
    call spamm_tree_print_leaves_2d_symm_recur(A%child_11, file_unit)

  end subroutine spamm_tree_print_leaves_2d_symm

END module spamm_elementals
