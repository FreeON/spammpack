!> Defines conversion operation between different data structures and SpAMM.
!!
!! @copyright
!!
!! Copyright (c) 2015, Los Alamos National Laboratory
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice, this
!! list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!! this list of conditions and the following disclaimer in the documentation
!! and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its contributors
!! may be used to endorse or promote products derived from this software without
!! specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! @author Matt Challacombe matt.challacombe@freeon.org
!! @author Nicolas Bock nicolasbock@freeon.org

module spamm_conversion

  use spamm_types
  use spamm_structors

  implicit none

contains

  FUNCTION SpAMM_convert_dense_to_tree_2d_symm( A, A_2d_O) RESUT(A_2d)

    real(SpAMM_KIND), dimension(:,:), intent(IN)    :: A
    type(SpAMM_tree_2d_symm) ,pointer,  optional    :: A_2d_O
    type(SpAMM_tree_2d_symm) ,pointer               :: A_2d
    
    IF(PRESENT(A_2d_O))THEN       
       a_2d => A_2d_O ! do this in place
    ELSE      
       ! lets get a fresh tree ... 
       a_2d => SpAMM_new_top_tree_2d_symm ( SIZE(A,1), SIZE(A,2) )
    ENDIF

    CALL SpAMM_convert_dense_to_tree_2d_symm_recur ( A, a_2d )

  END FUNCTION SpAMM_convert_dense_to_tree_2d_symm

  !> Recursively convert a dense matrix to a quadtree.
  RECURSIVE SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur ( A, A_2d )

    real(SpAMM_KIND), dimension(:,:), intent(IN)    :: A
    type(QuTree)    , pointer,        intent(INOUT) :: A_2d
    INTEGER, DIMENSION(2,2), pointer                :: bb

    M=SIZE(A,1)
    N=SIZE(A,2)

    bb => tree%frill%bndbx 

    ! peal off white space as its found 
    if(bb(0,1)>M)RETURN 
    if(bb(0,2)>N)RETURN 

    ! Leaf condition ? 
    wid=bb(1,1)-bb(0,1) ! w = [i]-[o]
    IF( wid == SPAMM_BLOCK_SIZE )THEN 

       if(.not.allocated(a_2d%chunk))allocate(a_2d%chunk(SPAMM_BLOCK_SIZE,SPAMM_BLOCK_SIZE))
       a_2d%chunk = SpAMM_Zero

       ! account for margin space ...
       m_marg = MIN( bb(1,1), M )
       n_marg = MIN( bb(1,2), N ) 

       ! data on the page ...    
       a_2d%chunk( 1:m_marg-bb(0,1)+1, 1:n_marg-bb(0,2)+1 ) = A( bb(0,1):m_marg, bb(0,2):n_marg ) 

    ELSE

       ! recur generically here ...
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_tree_2d_construct_00(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_tree_2d_construct_01(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_tree_2d_construct_11(A_2d) )

    ENDIF

    ! update the garnish 
    CALL SpAMM_redecorate_2d_symm(A_2d)

  END SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur

end module spamm_conversion



