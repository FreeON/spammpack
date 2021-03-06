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

  use spamm_algebra
  use spamm_globals
  use spamm_management
  use spamm_real_precision
  use spamm_types
  use spamm_utilities

#include "spamm_utility_macros.h"

  implicit none

contains

  !> Recursively convert a dense vector to a bitree.
  !!
  !! @param V The dense vector.
  !! @param qV A pointer to a bitree node.
  !! @param i_lower The lower value of the row index.
  !! @param i_upper The upper value of the row index.
  recursive subroutine spamm_convert_dense_to_bitree (V, qV, i_lower, i_upper)

    real(spamm_kind), dimension(:), intent(in) :: V
    type(bitree), pointer, intent(inout) :: qV
    integer, intent(in) :: i_lower, i_upper
    integer :: V_rows, convert_rows

    if(associated(qV)) then
       LOG_FATAL("qV should not already be associated")
       error stop
    endif

    LOG_DEBUG("converting: "//to_string(i_lower)//" "//to_string(i_upper))

    if(i_lower > size(V, 1)) then
       LOG_DEBUG("outside dense vector")
       return
    endif

    V_rows = i_upper-i_lower+1

    if(.not. associated(qV)) then
       LOG_DEBUG("allocating new node")
       call newbinode(qV, i_lower, i_upper)
    endif

    LOG_DEBUG("q: "//to_string(i_lower)//" "//to_string(i_upper))

    if(V_rows <= SPAMM_BLOCK_SIZE) then
       if(V_rows < SPAMM_BLOCK_SIZE) then
          LOG_FATAL("LOGIC ERROR IN SpAMM: padding error")
          LOG_FATAL("V_rows = "//to_string(V_rows))
          LOG_FATAL("SPAMM_BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE))
          error stop
       ELSE
          LOG_DEBUG("allocating new vect")

          if(allocated(qV%vect)) then
             deallocate(qV%vect)
          endif

          allocate(qV%vect(SPAMM_BLOCK_SIZE))

          ! Set new block to zero.
          qV%vect = 0

          ! Copy matrix elements.
          convert_rows = min(i_upper, size(V, 1))

          LOG_DEBUG("1:"//to_string(convert_rows-i_lower+1))
          LOG_DEBUG(to_string(i_lower)//":"//to_string(convert_rows))
          qV%vect(1:convert_rows-i_lower+1) = V(i_lower:convert_rows)
          LOG_DEBUG("stored: "//to_string(qV%vect))

          qV%norm = dot_product(qV%vect, qV%vect)

          ! Count number non-zeros.
          qV%number_nonzeros = count_nonzero(qV%vect)

          LOG_DEBUG("non-zeros: "//to_string(qV%number_nonzeros))
       endif
    else
       ! Avoid slicing here for performance.
       call spamm_convert_dense_to_bitree(V, qV%sect1, i_lower, i_lower+V_rows/2-1)
       call spamm_convert_dense_to_bitree(V, qV%sect2, i_lower+V_rows/2, i_upper)

       qV%number_nonzeros = 0
       qV%norm = 0

       if(associated(qV%sect1)) then
          qV%number_nonzeros = qV%number_nonzeros+qV%sect1%number_nonzeros
          qV%norm = qV%norm+qV%sect1%norm**2
       endif

       if(associated(qV%sect2)) then
          qV%number_nonzeros = qV%number_nonzeros+qV%sect2%number_nonzeros
          qV%norm = qV%norm+qV%sect2%norm**2
       endif

       qV%norm = sqrt(qV%norm)

    ENDIF

    LOG_DEBUG("done, going back up")

  end subroutine spamm_convert_dense_to_bitree

  !> Recursively convert a dense matrix to a quadtree.
  !!
  !! @param A The dense matrix.
  !! @param qA A pointer to a quadtree node.
  !! @param i_lower The lower value of the row index.
  !! @param i_upper The upper value of the row index.
  !! @param j_lower The lower value of the column index.
  !! @param j_upper The upper value of the column index.
  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree (A, qA, i_lower, i_upper, j_lower, j_upper)

    real(SpAMM_KIND), dimension(:,:), intent(IN) :: A
    type(QuTree), pointer, intent(INOUT) :: qA
    integer, intent(IN) :: i_lower, i_upper
    integer, intent(IN) :: j_lower, j_upper

    integer :: A_rows, A_cols
    integer :: convert_rows, convert_columns

    if(associated(qA)) then
       LOG_FATAL("qA should not already be associated")
       error stop
    endif

    LOG_DEBUG("converting: "//to_string(i_lower)//" "//to_string(i_upper))
    LOG_DEBUG("            "//to_string(j_lower)//" "//to_string(j_upper))

    if(i_lower > size(A, 1) .or. j_lower > size(A, 2)) then
       LOG_DEBUG("outside dense matrix")
       return
    endif

    A_rows = i_upper-i_lower+1
    A_cols = j_upper-j_lower+1

    if(.not. associated(qA)) then
       LOG_DEBUG("allocating new node")
       call newqunode(qA, i_lower, j_lower, i_upper, j_upper)
    endif

    LOG_DEBUG(to_string(qA))

    IF(A_rows <= SPAMM_BLOCK_SIZE .AND. A_cols <= SPAMM_BLOCK_SIZE)THEN
       IF(A_rows < SPAMM_BLOCK_SIZE .OR. A_cols < SPAMM_BLOCK_SIZE) THEN
          LOG_FATAL("[XgpSLv6M8u5ASgg3] LOGIC ERROR IN SpAMM: padding error")
          LOG_FATAL("A_rows = "//to_string(A_rows))
          LOG_FATAL("A_cols = "//to_string(A_cols))
          LOG_FATAL("SPAMM_BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE))
          error stop
       ELSE
          LOG_DEBUG("allocating new blok")

          if(allocated(qA%blok)) then
             deallocate(qA%blok)
          endif

          ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))

          ! Set new block to zero.
          qA%Blok = 0

          ! Copy matrix elements.
          convert_rows = min(i_upper, size(A, 1))
          convert_columns = min(j_upper, size(A, 2))

          qA%Blok(1:convert_rows-i_lower+1, 1:convert_columns-j_lower+1) = &
               A(i_lower:convert_rows, j_lower:convert_columns)

#ifdef SPAMM_STORE_TRANSPOSE
          if(allocated(qA%transpose_block)) then
             deallocate(qA%transpose_block)
          endif

          ALLOCATE(qA%transpose_block(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
          qA%transpose_block = transpose(qA%blok)
#endif

          qA%norm = sqrt(sum(qA%blok**2))

          ! Count number non-zeros.
          qA%number_nonzeros = count_nonzero(qA%blok)

          LOG_DEBUG("non-zeros: "//to_string(qA%number_nonzeros))
       ENDIF
    ELSE
       ! Avoid slicing here for performance.
       CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad11, &
            i_lower, &
            i_lower+A_rows/2-1, &
            j_lower, &
            j_lower+A_cols/2-1)
       CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad12, &
            i_lower, &
            i_lower+A_rows/2-1, &
            j_lower+A_cols/2, &
            j_upper)
       CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad21, &
            i_lower+A_rows/2, &
            i_upper, &
            j_lower, &
            j_lower+A_cols/2-1)
       CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad22, &
            i_lower+A_rows/2, &
            i_upper, &
            j_lower+A_cols/2, &
            j_upper)

       qA%number_nonzeros = 0
       qA%norm = 0

       if(associated(qA%quad11)) then
          qA%number_nonzeros = qA%number_nonzeros+qA%quad11%number_nonzeros
          qA%norm = qA%norm+qA%quad11%norm**2
       endif

       if(associated(qA%quad12)) then
          qA%number_nonzeros = qA%number_nonzeros+qA%quad12%number_nonzeros
          qA%norm = qA%norm+qA%quad12%norm**2
       endif

       if(associated(qA%quad21)) then
          qA%number_nonzeros = qA%number_nonzeros+qA%quad21%number_nonzeros
          qA%norm = qA%norm+qA%quad21%norm**2
       endif

       if(associated(qA%quad22)) then
          qA%number_nonzeros = qA%number_nonzeros+qA%quad22%number_nonzeros
          qA%norm = qA%norm+qA%quad22%norm**2
       endif

       qA%norm = sqrt(qA%norm)

    ENDIF

    LOG_DEBUG("done, going back up")

  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree

  !> Convert a dense 2nd order matrix to SpAMM matrix.
  !!
  !! @param A_dense The dense matrix.
  !!
  !! @result The SpAMM matrix.
  function spamm_convert_dense_to_matrix_2nd_order (A_dense) result (A)

    type(spamm_matrix_order_2), pointer :: A
    real(spamm_kind), dimension(:, :), intent(in) :: A_dense

    LOG_INFO("converting dense matrix")
    call spamm_allocate_matrix_2nd_order(size(A_dense, 1), size(A_dense, 2), A)
    call spamm_convert_dense_2_qutree(A_dense, A%root, 1, A%N_padded, 1, A%N_padded)

    if(associated(A%root)) then
       A%norm = A%root%norm
       A%number_nonzeros = A%root%number_nonzeros
    endif

    LOG_INFO("norm = "//to_string(A%norm))
    LOG_INFO("nnonzeros = "//to_string(A%number_nonzeros))
    LOG_INFO("done converting")

  end function spamm_convert_dense_to_matrix_2nd_order

  !> Convert a dense vector to a 1st order SpAMM object.
  !!
  !! @param V_dense The dense vector.
  !!
  !! @result The SpAMM vector.
  function spamm_convert_dense_to_order_1 (V_dense) result(V)

    real(spamm_kind), dimension(:), intent(in) :: V_dense
    type(spamm_matrix_order_1), pointer :: V

    LOG_INFO("converting dense vector")
    call spamm_allocate_matrix_order_1(size(V_dense, 1), V)
    call spamm_convert_dense_to_bitree(V_dense, V%root, 1, V%N_padded)

    if(associated(V%root)) then
       V%norm = V%root%norm
       V%number_nonzeros = V%root%number_nonzeros
    endif

    LOG_INFO("norm = "//to_string(V%norm))
    LOG_INFO("nnonzeros = "//to_string(V%number_nonzeros))
    LOG_INFO("done converting")

  end function spamm_convert_dense_to_order_1

  !> Convert a SpAMM matrix 2nd order to a dense matrix.
  !!
  !! @param A The SpAMM matrix.
  !! @param A_dense The dense matrix.
  subroutine spamm_convert_order_2_to_dense (A, A_dense)

    type(spamm_matrix_order_2), pointer, intent(in) :: A
    real(spamm_kind), dimension(:, :), allocatable, intent(inout) :: A_dense
    integer :: i, j

    LOG_DEBUG("printing matrix")

    if(allocated(A_dense)) then
       deallocate(A_dense)
    endif

    allocate(A_dense(A%M, A%N))

    do i = 1, A%M
       do j = 1, A%N
          A_dense(i, j) = get(A, i, j)
       enddo
    enddo

  end subroutine spamm_convert_order_2_to_dense

  !> Print a matrix.
  !!
  !! @param A The matrix.
  !! @param matrix_name The optional label.
  subroutine print_spamm_2nd_order (A, matrix_name)

    use spamm_utilities

    type(spamm_matrix_order_2), pointer, intent(in) :: A
    character(len = *), intent(in), optional :: matrix_name
    real(spamm_kind), dimension(:, :), allocatable :: A_dense
    character(len = 2000) :: format_string
    integer :: i, j

    call spamm_convert_order_2_to_dense(A, A_dense)

    write(format_string, "(A)") "("//to_string(size(A_dense, 2))//"ES11.3)"

    if(present(matrix_name)) then
       write(*, "(A)") "matrix "//trim(matrix_name)//" ="
    endif

    do i = 1, size(A_dense, 1)
       write(*, format_string) (A_dense(i, j), j = 1, size(A_dense, 2))
    enddo

  end subroutine print_spamm_2nd_order

  !> Recursively print a matrix tree.
  !!
  !! @param q A matrix node.
  recursive subroutine print_spamm_2nd_order_tree_node (q)

    type(qutree), pointer, intent(in) :: q

    if(.not. associated(q)) then
       return
    endif

    write(*, "(A)") to_string(q)

    if(associated(q%quad11)) then
       call print_spamm_2nd_order_tree_node(q%quad11)
    endif

    if(associated(q%quad12)) then
       call print_spamm_2nd_order_tree_node(q%quad12)
    endif

    if(associated(q%quad21)) then
       call print_spamm_2nd_order_tree_node(q%quad21)
    endif

    if(associated(q%quad22)) then
       call print_spamm_2nd_order_tree_node(q%quad22)
    endif

  end subroutine print_spamm_2nd_order_tree_node

  !> Print a matrix tree.
  !!
  !! @param A The matrix.
  !! @param matrix_name The optional label.
  subroutine print_spamm_2nd_order_tree (A, matrix_name)

    type(spamm_matrix_order_2), pointer, intent(in) :: A
    character(len = *), intent(in), optional :: matrix_name

    if(present(matrix_name)) then
       write(*, "(A)") "matrix tree "//trim(matrix_name)//" ="
    endif

    if(.not. associated(A)) then
       return
    endif

    if(.not. associated(A%root)) then
       write(*, "(A)") "root not associated"
       return
    endif

    call print_spamm_2nd_order_tree_node(A%root)

  end subroutine print_spamm_2nd_order_tree

end module spamm_convert
