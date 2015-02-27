module spamm_conversion

  use spamm_structures
  use spamm_xstructors

  implicit none

contains

  FUNCTION SpAMM_convert_dense_to_tree_2d_symm( A, A_2d_O) RESULT(A_2d)

    real(SpAMM_KIND), dimension(:,:), intent(IN)    :: A
    type(SpAMM_tree_2d_symm) ,pointer,  optional    :: A_2d_O
    type(SpAMM_tree_2d_symm) ,pointer               :: A_2d
    
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
    integer, dimension(0:1,1:2)                 :: bb
    INTEGER                                     :: M, N, M_marg, N_marg

    ! local variables ...
    M  =  a_2d%frill%NDimn(1)
    N  =  a_2d%frill%NDimn(2) 
    bb =  a_2d%frill%bndbx 

    ! peal off white space as found ... 
    if(bb(0,1)>M)RETURN 
    if(bb(0,2)>N)RETURN 

    if(a_2d%frill%leaf)then! Leaf condition ? 

       ! account for margin space ...
       M_marg = MIN( bb(1,1), M )
       N_marg = MIN( bb(1,2), N ) 
 
       ! move data on the page ...    
       a_2d%chunk( 1:m_marg-bb(0,1)+1, 1:n_marg-bb(0,2)+1 ) = A( bb(0,1):m_marg, bb(0,2):n_marg ) 

    ELSE ! recur generically here, poping with construct as needed ...
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_00(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_01(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_construct_tree_2d_symm_11(A_2d) )
    ENDIF

    ! update the garnish 
    CALL SpAMM_redecorate_tree_2d_symm(A_2d)
    write(*,*)' bb = ',a_2d%frill%bndbx

  END SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur

end module spamm_conversion



