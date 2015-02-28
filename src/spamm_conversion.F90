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

    WRITE(*,*)' a2d% NDIMN = ',a_2d%frill%ndimn


!    CALL SpAMM_convert_dense_to_tree_2d_symm_recur ( A, a_2d )
    WRITE(*,*)' a2d% NDIMN = ',a_2d%frill%ndimn

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
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_11(A_2d) )
    ENDIF

    ! update the garnish 
    CALL SpAMM_redecorate_tree_2d_symm(A_2d)
!    write(*,*)' bb = ',a_2d%frill%bndbx

  END SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur

end module spamm_conversion



