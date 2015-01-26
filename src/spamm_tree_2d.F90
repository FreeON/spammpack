!> Defines derived types used in SpAMMPACK.
!!
!! @copyright
!!
!! Copyright (c) 2014, Los Alamos National Laboratory
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

  !> The pointer to a tree node.
  type :: spamm_node_2d

     !> The actual tree node.
     type(tree_2d), pointer :: node => null()

  end type spamm_node_2d

  !> The SpAMM vector type.
  type :: tree_2d

     !> The tree node decoration.
     type(decoration_2d) :: decoration

     !> The Frobenius norm.
     real(SPAMM_KIND) :: norm = -1

     !> The number of non-zero elements.
     real(kind(0d0)) :: number_nonzeros = -1

     !> The axis-aligned bounding box.
     type(bounding_box_2d) :: bounding_box

     !> The pointers to the bisecting subtrees.
     type(spamm_node_2d) :: child(2)

     !> The vector data.
     real(SPAMM_KIND), dimension(:), allocatable :: data

   contains

#ifdef HAVE_FINALIZE
     !> The destructor.
     final :: delete_tree_2d
#endif

     !> Copy a vector.
     procedure :: copy_tree_2d_to_tree_2d

     !> The assignemnt operator.
     generic :: assignment(=) => copy_tree_2d_to_tree_2d

     !> String representation.
     procedure :: to_string => tree_2d_to_string

  end type tree_2d

#ifdef HAVE_CONSTRUCTOR
  !> The constructor interface.
  interface tree_2d
     module procedure new_tree_2d
     module procedure new_node_2d
  end interface tree_2d
#endif

  !> The destrcutor Interface.
  interface delete
     module procedure delete_tree_2d
  end interface delete

contains

  !> The constructor.
  !!
  !! @param N The matrix dimension.
  !!
  !! @return The tree node.
  type(tree_2d) function new_tree_2d (N) result(tree)

    use spamm_globals
    use spamm_utilities

    integer, intent(in) :: N

    integer :: i, N_padded, depth

    LOG_DEBUG("constructing new matrix")

    tree%norm = 0
    tree%number_nonzeros = 0

    i = 0
    do while(.true.)
       depth = i
       N_padded = SPAMM_BLOCK_SIZE*2**i
       if(N_padded > N) then
          exit
       endif
       i = i+1
    enddo

    tree%decoration = decoration_2d(N, N_padded, depth)

#ifdef HAVE_CONSTRUCTOR
    tree%bounding_box = bounding_box_2d(N_padded/2, N_padded/2)
#else
    tree%bounding_box = new_bounding_box_2d(N_padded/2, N_padded/2)
#endif

  end function new_tree_2d

  !> The constructor.
  !!
  !! @param N The matrix dimension.
  !! @param N_padded The padded matrix dimension.
  !! @param depth The tree depth.
  !! @param bounding_box The axis-aligned bounding box.
  !!
  !! @return The tree node.
  type(tree_2d) function new_node_2d (decoration, bounding_box) result(tree)

    use spamm_globals
    use spamm_utilities

    type(decoration_2d), intent(in) :: decoration
    type(bounding_box_2d), intent(in) :: bounding_box

    LOG_DEBUG("constructing new matrix")

    tree%decoration = decoration
    tree%bounding_box = bounding_box

    tree%norm = 0
    tree%number_nonzeros = 0

  end function new_node_2d

  !> The destructor.
  !!
  !! @param self The object to deallocate.
  elemental subroutine delete_tree_2d (self)

    type(tree_2d), intent(inout) :: self

    if(allocated(self%data)) then
       deallocate(self%data)
    endif

  end subroutine delete_tree_2d

  !> Copy a vector to another vector:
  !!
  !! \f$ A \leftarrow B \f$.
  !!
  !! @param A The vector A.
  !! @param B The vector B.
  recursive subroutine copy_tree_2d_to_tree_2d (A, B)

    use spamm_utilities

    class(tree_2d), intent(out) :: A
    class(tree_2d), intent(in) :: B

    integer :: i

    LOG_DEBUG("copying vector")

    A%decoration = B%decoration
    A%norm = B%norm
    A%number_nonzeros = B%number_nonzeros
    A%bounding_box = B%bounding_box

    LOG_DEBUG("here...")

    if(allocated(B%data)) then
       LOG_DEBUG("copying matrix data")
       A%data = B%data
    else
       do i = 1, 2
          if(associated(B%child(i)%node)) then
             LOG_DEBUG("descending")
             call copy_tree_2d_to_tree_2d(A%child(i)%node, B%child(i)%node)
          endif
       enddo
    endif

    LOG_DEBUG("done copying")

  end subroutine copy_tree_2d_to_tree_2d

  !> String representation of matrix.
  !!
  !! @param A The matrix.
  !!
  !! @return The string representation.
  character(len = 1000) function tree_2d_to_string (A) result(string)

    class(tree_2d), intent(in) :: A

    character(len = 100) :: temp

    write(string, "(A)") "N = "//trim(adjustl(A%decoration%to_string()))

    write(temp, "(ES15.5)") A%norm
    write(string, "(A)") trim(string)//", norm = "//trim(adjustl(temp))

    write(temp, "(ES15.5)") A%number_nonzeros
    write(string, "(A)") trim(string)//", nnonz = "//trim(adjustl(temp))

    write(temp, *) A%bounding_box%left_edge()
    write(string, "(A)") trim(string)//", bbox = [ "//trim(adjustl(temp))

    write(temp, *) A%bounding_box%right_edge()
    write(string, "(A)") trim(string)//", "//trim(adjustl(temp))//" ]"

    if(allocated(A%data)) then
       write(string, "(A)") trim(string)//", data is allocated"
    endif

  end function tree_2d_to_string

end module spamm_tree_2d
