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
  use spamm_real_precision

  implicit none

  !> Symmetric matrix.
  type :: tree_2d_symmetric

     !> The tree node decoration.
     type(decoration_2d) :: decoration

     !> The axis-aligned bounding box.
     type(bounding_box_2d) :: bounding_box

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

  !> General matrix.
  type :: tree_2d_general

     !> The tree node decoration.
     type(decoration_2d) :: decoration

     !> The axis-aligned bounding box.
     type(bounding_box_2d) :: bounding_box

     !> The vector data.
     real(SPAMM_KIND), allocatable :: data(:, :)

     !> Pointer to upper left quadrant.
     type(tree_2d_general), pointer :: child_00 => null()

     !> Pointer to upper left quadrant.
     type(tree_2d_general), pointer :: child_01 => null()

     !> Pointer to upper left quadrant.
     type(tree_2d_general), pointer :: child_10 => null()

     !> Pointer to upper left quadrant.
     type(tree_2d_general), pointer :: child_11 => null()

  end type tree_2d_general

contains

  !> The constructor.
  !!
  !! @param N The matrix dimension.
  !!
  !! @return The tree node.
  function new_tree_2d_symmetric (N) result (tree)

    use spamm_globals
    use spamm_strings

    type(tree_2d_symmetric), pointer :: tree
    integer, intent(in) :: N

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

    tree%decoration = decoration_2d(N, N_padded, depth, 0, 0)
    tree%bounding_box = bounding_box_2d(1, N_padded, 1, N_padded)

  end function new_tree_2d_symmetric

  function new_tree_2d_symmetric_decorated (decoration) result (tree)

    type(tree_2d_symmetric), pointer :: tree
    type(decoration_2d), intent(in) :: decoration

    tree => new_tree_2d_symmetric(decoration%N)

  end function new_tree_2d_symmetric_decorated

  !> The identity matrix.
  !!
  !! @param N The matrix dimension.
  !! @return The new matrix.
  function identity_tree_2d_symmetric (N) result (tree)

    type(tree_2d_symmetric), pointer :: tree
    integer, intent(in) :: N

    tree => new_tree_2d_symmetric(N)
    call set_identity_2d_symmetric(tree)

  end function identity_tree_2d_symmetric

  recursive subroutine set_identity_2d_symmetric (tree)

    use spamm_bisect
    use spamm_globals

    type(tree_2d_symmetric), intent(inout) :: tree

    integer :: row_middle, column_middle
    integer :: i

    row_middle = bisect(tree%bounding_box%row_lower, tree%bounding_box%row_upper)
    column_middle = bisect(tree%bounding_box%column_lower, tree%bounding_box%column_upper)

    if(row_middle < 0 .or. column_middle < 0) then
       ! Leaf node.
       allocate(tree%data(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
       do i = 1, SPAMM_BLOCK_SIZE
          tree%data(i, i) = 1
       enddo
    else
       ! Recur down the diagonal.
       tree%child_00 => new_tree_2d_symmetric_decorated(tree%decoration)
       tree%child_11 => new_tree_2d_symmetric_decorated(tree%decoration)

       tree%child_00%bounding_box = bounding_box_2d(tree%bounding_box%row_lower, &
            row_middle, &
            tree%bounding_box%column_lower, &
            column_middle)

       tree%child_11%bounding_box = bounding_box_2d(row_middle+1, &
            tree%bounding_box%row_upper, &
            column_middle+1, &
            tree%bounding_box%column_upper)

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
