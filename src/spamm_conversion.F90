module spamm_conversion

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals

  implicit none

  interface spamm_convert_from_dense
     module procedure spamm_convert_dense_to_tree_2d_symm
  end interface spamm_convert_from_dense

contains

  !> Convert a dense matrix to a quadtree.
  !!
  !! @param a The dense matrix.
  !! @param a_old An already existing tree. Will be updated in place.
  !! @return The converted matrix as a quadtree.
  function spamm_convert_dense_to_tree_2d_symm(a, a_old) result(a_2d)

    real(SPAMM_KIND), dimension(:,:), intent(in) :: a
    type(spamm_tree_2d_symm), pointer, optional :: a_old
    type(spamm_tree_2d_symm), pointer :: a_2d

    if(present(a_old)) then
       a_2d => a_old ! Data pass in, keep it in place.
    end if

    if(.not.associated(a_2d)) then
       a_2d => spamm_new_top_tree_2d_symm([size(a,1), size(a,2)]) ! A new tree.
    end if

    call spamm_flip(a_2d)
    call spamm_convert_dense_to_tree_2d_symm_recur(a, a_2d)
    call spamm_prune(a_2d)

  end function spamm_convert_dense_to_tree_2d_symm

  !> Recursively convert a dense matrix to a quadtree.
  recursive subroutine spamm_convert_dense_to_tree_2d_symm_recur(a, a_2d)

    real(SPAMM_KIND), dimension(:,:), intent(in) :: a
    type(spamm_tree_2d_symm), pointer :: a_2d

    integer, dimension(1:2) :: lo, hi

    if(.not.associated(a_2d)) return

    if(a_2d%frill%leaf) then! Leaf condition ?
       a_2d%frill%needs_initialization = .false.
       lo = a_2d%frill%bndbx(0,:)
       hi = a_2d%frill%bndbx(1,:)
       ! move data on the page ...
       a_2d%chunk(1:(hi(1)-lo(1)+1), 1:(hi(2)-lo(2)+1)) = A(lo(1):hi(1), lo(2):hi(2))
    else ! recur generically here, coping with construct as needed ...
       call spamm_convert_dense_to_tree_2d_symm_recur(a, spamm_construct_tree_2d_symm_00(a_2d))
       call spamm_convert_dense_to_tree_2d_symm_recur(a, spamm_construct_tree_2d_symm_01(a_2d))
       call spamm_convert_dense_to_tree_2d_symm_recur(a, spamm_construct_tree_2d_symm_10(a_2d))
       call spamm_convert_dense_to_tree_2d_symm_recur(a, spamm_construct_tree_2d_symm_11(a_2d))
    end if

    ! update the garnish
    call spamm_redecorate_tree_2d_symm(a_2d)

  end subroutine spamm_convert_dense_to_tree_2d_symm_recur

  subroutine SpAMM_convert_tree_2d_symm_to_dense(A_2d, A)

    type(SpAMM_tree_2d_symm),         pointer           :: A_2d
    real(SPAMM_KIND), dimension(:,:)                    :: A

    !    IF(.not.allocated(A))THEN
    !       STOP' need pre allocation here ...'
    !       ALLOCATE(A( 1:A_2d%frill%bndbx(1,1),  &
    !            1:A_2d%frill%bndbx(1,2) ))
    !   END IF

    A=0
    call SpAMM_convert_tree_2d_symm_to_dense_recur (A_2d,A)

  end subroutine SpAMM_convert_tree_2d_symm_to_dense

  !> Recursively convert a dense matrix to a quadtree.
  recursive subroutine SpAMM_convert_tree_2d_symm_to_dense_recur (A_2d,A)

    real(SPAMM_KIND), dimension(:,:)     :: A
    !    real(SPAMM_KIND), dimension(:,:)     :: A
    type(SpAMM_tree_2d_symm),        pointer :: A_2d
    integer, dimension(1:2)                  :: lo,hi

    if(.not.associated(a_2d))return

    if(a_2d%frill%leaf)then! Leaf condition ?

       lo=a_2d%frill%bndbx(0,:)
       hi=a_2d%frill%bndbx(1,:)

       ! move data on the page ...
       A(lo(1):hi(1),lo(2):hi(2))=a_2d%chunk( 1:(hi(1)-lo(1)+1), 1:(hi(2)-lo(2)+1))

    else ! recur generically here ...
       call SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_00, A )
       call SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_01, A )
       call SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_10, A )
       call SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_11, A )
    end if
  end subroutine SpAMM_convert_tree_2d_symm_to_dense_recur

end module spamm_conversion
