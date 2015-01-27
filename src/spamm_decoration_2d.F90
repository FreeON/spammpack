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
module spamm_decoration_2d

#include "spamm_utility_macros.h"

  implicit none

  type :: decoration_2d

     !> The number of entries.
     integer :: N = -1

     !> The padded vector dimension.
     integer :: N_padded = -1

     !> The tree depth. The root tier is 0, tier == depth is the leaf
     !> node tier, i.e. the tier at which actual matrix elements are
     !> stored.
     integer :: depth = -1

   contains

     procedure :: to_string => decoration_2d_to_string

  end type decoration_2d

contains

  !> String representation of decoration.
  !!
  !! @param decoration The node decoration.
  !!
  !! @return The string representation.
  character(len = 1000) function decoration_2d_to_string (decoration) result(string)

    class(decoration_2d), intent(in) :: decoration

    character(len = 100) :: temp

    write(temp, *) decoration%N
    write(string, "(A)") "N = "//trim(adjustl(temp))

    write(temp, *) decoration%N_padded
    write(string, "(A)") trim(string)//", N_padded = "//trim(adjustl(temp))

    write(temp, *) decoration%depth
    write(string, "(A)") trim(string)//", depth = "//trim(adjustl(temp))

  end function decoration_2d_to_string

end module spamm_decoration_2d
