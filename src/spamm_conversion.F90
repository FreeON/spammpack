module spamm_conversion

  use spamm_structures
  use spamm_xstructors

  implicit none

contains

  FUNCTION SpAMM_convert_dense_to_tree_2d_symm( A, A_2d_O) RESULT(A_2d)

    real(SpAMM_KIND), dimension(:,:), intent(IN)    :: A
    type(SpAMM_tree_2d_symm) ,pointer,  optional    :: A_2d_O
    type(SpAMM_tree_2d_symm) ,pointer               :: A_2d
    
    WRITE(*,*)' SIZE A = ',(/ SIZE(A,1), SIZE(A,2) /)

    IF(PRESENT(A_2d_O))THEN       
       a_2d => A_2d_O ! do this in place
    ELSE      
       ! lets get a fresh tree ... 
       a_2d => SpAMM_new_top_tree_2d_symm ((/ SIZE(A,1), SIZE(A,2) /))
    ENDIF

    CALL SpAMM_convert_dense_to_tree_2d_symm_recur ( A, a_2d )

  END FUNCTION SpAMM_convert_dense_to_tree_2d_symm

  !> Recursively convert a dense matrix to a quadtree.
  RECURSIVE SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur (A,A_2d)

    real(SpAMM_KIND), dimension(:,:),intent(in) :: A
    type(SpAMM_tree_2d_symm), pointer           :: A_2d
    integer, dimension(1:2)                     :: lo,hi

    if(.not.associated(a_2d))return
    if(a_2d%frill%leaf)then! Leaf condition ? 

       lo=a_2d%frill%bndbx(0,:) 
       hi=a_2d%frill%bndbx(1,:) 
 
       ! move data on the page ...    
       a_2d%chunk( 1:(hi(1)-lo(1)+1), 1:(hi(2)-lo(2)+1))=A(lo(1):hi(1),lo(2):hi(2))

    ELSE ! recur generically here, poping with construct as needed ...

       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_00(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_01(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_10(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_11(A_2d) )

    ENDIF

    ! update the garnish 
    CALL SpAMM_redecorate_tree_2d_symm(A_2d)

  END SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur




  FUNCTION SpAMM_convert_tree_2d_symm_to_dense(A_2d, A_O) RESULT(A)

    type(SpAMM_tree_2d_symm),         pointer           :: A_2d
    real(SpAMM_KIND), dimension(:,:), pointer, optional :: A_O
    real(SpAMM_KIND), dimension(:,:), pointer           :: A

    IF(PRESENT(A_O))THEN       
       A=>A_O 
    ELSE      
       ALLOCATE(A( 1:A_2d%frill%bndbx(1,1),  &
                   1:A_2d%frill%bndbx(1,2) ))
    ENDIF

    A=SpAMM_zero
    CALL SpAMM_convert_tree_2d_symm_to_dense_recur (A_2d,A)

  END FUNCTION SpAMM_convert_tree_2d_symm_to_dense

  !> Recursively convert a dense matrix to a quadtree.
  RECURSIVE SUBROUTINE SpAMM_convert_tree_2d_symm_to_dense_recur (A_2d,A)

    real(SpAMM_KIND), dimension(:,:),pointer :: A
    type(SpAMM_tree_2d_symm),        pointer :: A_2d
    integer, dimension(1:2)                  :: lo,hi

    if(.not.associated(a_2d))return

    if(a_2d%frill%leaf)then! Leaf condition ? 

       lo=a_2d%frill%bndbx(0,:) 
       hi=a_2d%frill%bndbx(1,:)  

       ! move data on the page ...    
       A(lo(1):hi(1),lo(2):hi(2))=a_2d%chunk( 1:(hi(1)-lo(1)+1), 1:(hi(2)-lo(2)+1))

    ELSE ! recur generically here ...
       CALL SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_00, A )
       CALL SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_01, A )
       CALL SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_10, A )
       CALL SpAMM_convert_tree_2d_symm_to_dense_recur( A_2d%child_11, A )
    ENDIF
  END SUBROUTINE SpAMM_convert_tree_2d_symm_to_dense_recur

end module spamm_conversion



