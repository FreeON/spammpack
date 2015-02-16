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

module spamm_convert

  use spamm_types
  use spamm_structors

  implicit none

contains

  FUNCTION SpAMM_convert_dense_to_tree_2d_symm( A, A_2d_O) RESUT(A_2d)

    real(SpAMM_KIND), dimension(:,:), intent(IN)    :: A
    type(SpAMM_tree_2d_symm) ,pointer,  optional    :: A_2d_O
    type(SpAMM_tree_2d_symm) ,pointer               :: A_2d
    
    IF(PRESENT(A_2d_O))THEN
       ! for whatever reason, lets do this in place
       a_2d => A_2d_O 
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

    ! peal off the white space as its found 
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
       a_2d%chunk( bb(0,1):m_marg, bb(0,2):n_marg ) = A( bb(0,1):m_marg, bb(0,2):n_marg ) 

    ELSE

       ! recur generically here ...
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_tree_2d_construct_00(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_tree_2d_construct_01(A_2d) )
       CALL SpAMM_convert_dense_to_tree_2d_symm_recur( A, SpAMM_tree_2d_construct_11(A_2d) )

    ENDIF

    CALL SpAMM_redecorate_2d_symm(A_2d)

  END SUBROUTINE SpAMM_convert_dense_to_tree_2d_symm_recur

end module spamm_convert


!!$  !> Recursively convert a dense vector to a bitree.
!!$  !!
!!$  !! @param V The dense vector.
!!$  !! @param qV A pointer to a bitree node.
!!$  !! @param i_lower The lower value of the row index.
!!$  !! @param i_upper The upper value of the row index.
!!$  recursive subroutine spamm_convert_dense_to_bitree (V, qV, i_lower, i_upper)
!!$
!!$    real(spamm_kind), dimension(:), intent(in) :: V
!!$    type(bitree), pointer, intent(inout) :: qV
!!$    integer, intent(in) :: i_lower, i_upper
!!$    integer :: V_rows, convert_rows
!!$
!!$    if(associated(qV)) then
!!$       LOG_FATAL("qV should not already be associated")
!!$       error stop
!!$    endif
!!$
!!$    LOG_DEBUG("converting: "//to_string(i_lower)//" "//to_string(i_upper))
!!$
!!$    if(i_lower > size(V, 1)) then
!!$       LOG_DEBUG("outside dense vector")
!!$       return
!!$    endif
!!$
!!$    V_rows = i_upper-i_lower+1
!!$
!!$    if(.not. associated(qV)) then
!!$       LOG_DEBUG("allocating new node")
!!$       call newbinode(qV, i_lower, i_upper)
!!$    endif
!!$
!!$    LOG_DEBUG("q: "//to_string(i_lower)//" "//to_string(i_upper))
!!$
!!$    if(V_rows <= SPAMM_BLOCK_SIZE) then
!!$       if(V_rows < SPAMM_BLOCK_SIZE) then
!!$          LOG_FATAL("LOGIC ERROR IN SpAMM: padding error")
!!$          LOG_FATAL("V_rows = "//to_string(V_rows))
!!$          LOG_FATAL("SPAMM_BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE))
!!$          error stop
!!$       ELSE
!!$          LOG_DEBUG("allocating new vect")
!!$
!!$          if(allocated(qV%vect)) then
!!$             deallocate(qV%vect)
!!$          endif
!!$
!!$          allocate(qV%vect(SPAMM_BLOCK_SIZE))
!!$
!!$          ! Set new block to zero.
!!$          qV%vect = 0
!!$
!!$          ! Copy matrix elements.
!!$          convert_rows = min(i_upper, size(V, 1))
!!$
!!$          LOG_DEBUG("1:"//to_string(convert_rows-i_lower+1))
!!$          LOG_DEBUG(to_string(i_lower)//":"//to_string(convert_rows))
!!$          qV%vect(1:convert_rows-i_lower+1) = V(i_lower:convert_rows)
!!$          LOG_DEBUG("stored: "//to_string(qV%vect))
!!$
!!$          qV%norm = dot_product(qV%vect, qV%vect)
!!$
!!$          ! Count number non-zeros.
!!$          qV%number_nonzeros = count_nonzero(qV%vect)
!!$
!!$          LOG_DEBUG("non-zeros: "//to_string(qV%number_nonzeros))
!!$       endif
!!$    else
!!$       ! Avoid slicing here for performance.
!!$       call spamm_convert_dense_to_bitree(V, qV%sect1, i_lower, i_lower+V_rows/2-1)
!!$       call spamm_convert_dense_to_bitree(V, qV%sect2, i_lower+V_rows/2, i_upper)
!!$
!!$       qV%number_nonzeros = 0
!!$       qV%norm = 0
!!$
!!$       if(associated(qV%sect1)) then
!!$          qV%number_nonzeros = qV%number_nonzeros+qV%sect1%number_nonzeros
!!$          qV%norm = qV%norm+qV%sect1%norm**2
!!$       endif
!!$
!!$       if(associated(qV%sect2)) then
!!$          qV%number_nonzeros = qV%number_nonzeros+qV%sect2%number_nonzeros
!!$          qV%norm = qV%norm+qV%sect2%norm**2
!!$       endif
!!$
!!$       qV%norm = sqrt(qV%norm)
!!$
!!$    ENDIF
!!$
!!$    LOG_DEBUG("done, going back up")
!!$
!!$  end subroutine spamm_convert_dense_to_bitree


end module spamm_convert
