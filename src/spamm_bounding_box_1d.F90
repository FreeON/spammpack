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
module spamm_bounding_box_1d

#include "spamm_utility_macros.h"

  use spamm_real_precision
  use spamm_utilities

  implicit none

  !> A bounding box.
  type :: bounding_box_1d

     !> The center.
     integer :: center

     !> The width.
     integer :: width

   contains

     !> Get the left edge.
     procedure :: left_edge

     !> Get the right edge.
     procedure :: right_edge

  end type bounding_box_1d

  !> The constructor.
  !interface new
  !   module procedure new_bounding_box_1d
  !end interface new

contains

  !> The constructor.
  function new_bounding_box_1d (center, width) result(box)

    type(bounding_box_1d), pointer :: box
    integer, intent(in) :: center, width

    LOG_DEBUG("constructing new bounding box")

    allocate(box)

    box%center = center
    box%width = width

  end function new_bounding_box_1d

  !> Get the left edge.
  !!
  !! @return The left edge of the bounding box.
  integer function left_edge (self)

    class(bounding_box_1d), intent(in) :: self
    left_edge = self%center-self%width+1

  end function left_edge

  !> Get the right edge.
  !!
  !! @return The right edge of the bounding box.
  integer function right_edge (self)

    class(bounding_box_1d), intent(in) :: self
    right_edge = self%center+self%width

  end function right_edge

end module spamm_bounding_box_1d
