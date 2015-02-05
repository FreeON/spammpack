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
     real(SPAMM_KIND), dimension(:), allocatable :: data

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
     real(SPAMM_KIND), dimension(:), allocatable :: data

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
  function new_tree_2d_symmetric (N) result(tree)

    use spamm_globals
    use spamm_utilities

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

    character(len = 100) :: temp

    write(string, "(A)") "N = "//trim(adjustl(A%decoration%to_string()))

    if(allocated(A%data)) then
       write(string, "(A)") trim(string)//", data is allocated"
    endif

  end function tree_2d_symmetric_to_string

end module spamm_tree_2d
