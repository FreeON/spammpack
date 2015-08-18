!> Defines derived types used in SpAMMPACK.
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
module spamm_tree_2d

#include "spamm_utility_macros.h"

  use spamm_bounding_box_2d
  use spamm_decoration_2d
  use spamm_extra_2d
  use spamm_real_precision

  implicit none

  !> Symmetric matrix.
  type :: tree_2d_symmetric

     !> The tree node decoration.
     type(decoration_2d) :: decoration

     !> Some extra information.
     type(extra_2d) :: extra

     !> The vector data.
     real(SPAMM_KIND), allocatable :: data(:, :)

     !> Pointer to upper left quadrant.
     type(tree_2d_symmetric), pointer :: child_00 => null()

     !> Pointer to upper left quadrant.
     type(tree_2d_symmetric), pointer :: child_01 => null()

     !> Pointer to upper left quadrant.
     type(tree_2d_symmetric), pointer :: child_11 => null()

   contains

     !> String representation.
     procedure :: to_string => tree_2d_symmetric_to_string

  end type tree_2d_symmetric

contains

  !> The constructor.
  !!
  !! @param M The number of rows.
  !! @param N The number of columns.
  !!
  !! @return The tree node.
  function new_tree_2d_symmetric (M, N) result (tree)

    use spamm_globals
    use spamm_strings

    type(tree_2d_symmetric), pointer :: tree
    integer, intent(in) :: M, N

    integer :: i, N_padded, depth

    LOG_DEBUG("constructing new symmetric matrix")

    allocate(tree)

    i = 0
    do while(.true.)
       depth = i
       N_padded = SPAMM_BLOCK_SIZE*2**i
       if(N_padded > N) then
          exit
       endif
       i = i+1
    enddo

    tree%decoration = decoration_2d(0, 0, bounding_box_2d([1, 1], [N_padded, N_padded]))
    tree%extra = extra_2d(M, N)

  end function new_tree_2d_symmetric

  function new_node_2d_symmetric (decoration, extra) result (tree)

    type(tree_2d_symmetric), pointer :: tree
    type(decoration_2d), intent(in) :: decoration
    type(extra_2d), intent(in) :: extra

    allocate(tree)

    tree%decoration = decoration
    tree%extra = extra

  end function new_node_2d_symmetric

  !> The identity matrix.
  !!
  !! @param N The matrix dimension.
  !! @return The new matrix.
  function identity_tree_2d_symmetric (N) result (tree)

    type(tree_2d_symmetric), pointer :: tree
    integer, intent(in) :: N

    tree => new_tree_2d_symmetric(N, N)
    call set_identity_2d_symmetric(tree)

  end function identity_tree_2d_symmetric

  recursive subroutine set_identity_2d_symmetric (tree)

    use spamm_bisect
    use spamm_globals

    type(tree_2d_symmetric), intent(inout) :: tree

    integer :: middle(0:1)
    integer :: i

    LOG_DEBUG("[identity] "//trim(tree%decoration%to_string()))

    middle(0) = bisect(tree%decoration%bounding_box%lower(0), tree%decoration%bounding_box%upper(0))
    middle(1) = bisect(tree%decoration%bounding_box%lower(1), tree%decoration%bounding_box%upper(1))

    if(middle(0) < 0 .or.middle(1) < 0) then
       ! Leaf node.
       allocate(tree%data(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
       do i = 1, SPAMM_BLOCK_SIZE
          tree%data(i, i) = 1
       enddo
    else
       ! Recur down the diagonal.
       tree%child_00 => new_node_2d_symmetric(decoration_2d(0, 0, &
            bounding_box_2d( &
            [tree%decoration%bounding_box%lower(0), tree%decoration%bounding_box%lower(1)], &
            [middle(0), middle(1)])), tree%extra)

       tree%child_11 => new_node_2d_symmetric(decoration_2d(0, 0, &
            bounding_box_2d([middle(0)+1, middle(1)+1], &
            [tree%decoration%bounding_box%upper(0), tree%decoration%bounding_box%upper(1)])), tree%extra)

       call set_identity_2d_symmetric(tree%child_00)
       call set_identity_2d_symmetric(tree%child_11)
    endif

  end subroutine set_identity_2d_symmetric

  !> The destructor.
  !!
  !! @param self The object to deallocate.
  recursive subroutine delete_tree_2d_symmetric (self)

    type(tree_2d_symmetric), pointer, intent(inout) :: self

    if(allocated(self%data)) then
       deallocate(self%data)
    endif

    deallocate(self)
    nullify(self)

  end subroutine delete_tree_2d_symmetric

  !> String representation of matrix.
  !!
  !! @param A The matrix.
  !!
  !! @return The string representation.
  function tree_2d_symmetric_to_string (A) result(string)

    character(len = 1000) :: string
    class(tree_2d_symmetric), intent(in) :: A

    write(string, "(A)") "N = "//trim(adjustl(A%decoration%to_string()))

    if(allocated(A%data)) then
       write(string, "(A)") trim(string)//", data is allocated"
    endif

  end function tree_2d_symmetric_to_string

end module spamm_tree_2d
