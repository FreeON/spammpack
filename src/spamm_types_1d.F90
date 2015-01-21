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
module spamm_types_1d

#include "spamm_utility_macros.h"

  use spamm_real_precision

  implicit none

  !> The pointer to a tree node.
  type :: spamm_node_1d

     !> The actual tree node.
     type(spamm_matrix_1d), pointer :: node => null()

  end type spamm_node_1d

  !> A bounding box.
  type :: bounding_box_1d

     !> The center.
     integer :: center

     !> The width.
     integer :: width

  end type bounding_box_1d

  !> The SpAMM vector type.
  type :: spamm_matrix_1d

     !> The number of entries.
     integer :: N = -1

     !> The padded vector dimension.
     integer :: N_padded = -1

     !> The tree depth. The root tier is 0, tier == depth is the leaf
     !> node tier, i.e. the tier at which actual matrix elements are
     !> stored.
     integer :: depth = -1

     !> The Frobenius norm.
     real(SPAMM_KIND) :: norm = -1

     !> The number of non-zero elements.
     real(kind(0d0)) :: number_nonzeros = -1

     !> The axis-aligned bounding box.
     type(bounding_box_1d) :: bounding_box

     !> The pointers to the bisecting subtrees.
     type(spamm_node_1d) :: child(2)

     !> The vector data.
     real(SPAMM_KIND), dimension(:), allocatable :: data

   contains

#ifdef HAVE_FINALIZE
     !> The destructor.
     final :: delete_matrix_1d
#endif

     !> Copy a vector.
     procedure :: copy_matrix_1d_to_matrix_1d

     !> The assignemnt operator.
     generic :: assignment(=) => copy_matrix_1d_to_matrix_1d

     !> String representation.
     procedure :: to_string

  end type spamm_matrix_1d

#ifdef HAVE_CONSTRUCTOR
  !> The constructor interface.
  interface spamm_matrix_1d
     module procedure new_matrix_1d
  end interface spamm_matrix_1d
#endif

  !> The destrcutor Interface.
  interface delete
     module procedure delete_matrix_1d
  end interface delete

#ifdef HAVE_CONSTRUCTOR
  !> The constructor.
  interface bounding_box_1d
     module procedure new_bounding_box_1d
  end interface bounding_box_1d
#endif

contains

  !> The constructor.
  type(bounding_box_1d) function new_bounding_box_1d (center, width) result(box)

    integer, intent(in) :: center, width

    box%center = center
    box%width = width

  end function new_bounding_box_1d

  !> The constructor.
  !!
  !! @param N The matrix dimension.
  type(spamm_matrix_1d) function new_matrix_1d (N) result(tree)

    use spamm_globals

    integer, intent(in) :: N

    integer :: i

    LOG_DEBUG("constructing new matrix")

    tree%N = N
    tree%norm = 0
    tree%number_nonzeros = 0

    i = 0
    do while(.true.)
       tree%depth = i
       tree%N_padded = SPAMM_BLOCK_SIZE*2**i
       if(tree%N_padded > N) then
          exit
       endif
       i = i+1
    enddo

    tree%bounding_box = bounding_box_1d(tree%N_padded/2, tree%N_padded/2)

  end function new_matrix_1d

  !> The destructor.
  !!
  !! @param self The object to deallocate.
  elemental subroutine delete_matrix_1d (self)

    type(spamm_matrix_1d), intent(inout) :: self

    if(allocated(self%data)) then
       deallocate(self%data)
    endif

  end subroutine delete_matrix_1d

  !> Copy a vector to another vector:
  !!
  !! \f$ A \leftarrow B \f$.
  !!
  !! @param A The vector A.
  !! @param B The vector B.
  recursive subroutine copy_matrix_1d_to_matrix_1d (A, B)

    class(spamm_matrix_1d), intent(out) :: A
    class(spamm_matrix_1d), intent(in) :: B

    integer :: i

    LOG_DEBUG("copying vector")

    A%N = B%N
    A%N_padded = B%N_padded
    A%depth = B%depth
    A%norm = B%norm
    A%number_nonzeros = B%number_nonzeros
    A%bounding_box = B%bounding_box

    if(allocated(B%data)) then
       LOG_DEBUG("copying matrix data")
       A%data = B%data
    else
       do i = 1, 2
          if(associated(B%child(i)%node)) then
             LOG_DEBUG("descending")
             call copy_matrix_1d_to_matrix_1d(A%child(i)%node, B%child(i)%node)
          endif
       enddo
    endif

  end subroutine copy_matrix_1d_to_matrix_1d

  !> String representation of matrix.
  !!
  !! @param A The matrix.
  !!
  !! @return The string representation.
  character(len = 1000) function to_string (A)

    class(spamm_matrix_1d), intent(in) :: A

    character(len = 100) :: temp

    write(temp, *) A%N
    write(to_string, "(A)") "N = "//trim(adjustl(temp))

    write(temp, *) A%N_padded
    write(to_string, "(A)") trim(to_string)//", N_padded = "//trim(adjustl(temp))

    write(temp, *) A%depth
    write(to_string, "(A)") trim(to_string)//", depth = "//trim(adjustl(temp))

    write(temp, "(ES20.10)") A%norm
    write(to_string, "(A)") trim(to_string)//", norm = "//trim(adjustl(temp))

    write(temp, "(ES20.10)") A%number_nonzeros
    write(to_string, "(A)") trim(to_string)//", nnonz = "//trim(adjustl(temp))

    write(temp, *) A%bounding_box%center-A%bounding_box%width+1
    write(to_string, "(A)") trim(to_string)//", bbox = [ "//trim(adjustl(temp))

    write(temp, *) A%bounding_box%center+A%bounding_box%width-1
    write(to_string, "(A)") trim(to_string)//", "//trim(adjustl(temp))//" ]"

    if(allocated(A%data)) then
       write(to_string, "(A)") trim(to_string)//", data is allocated"
    endif

  end function to_string

end module spamm_types_1d
