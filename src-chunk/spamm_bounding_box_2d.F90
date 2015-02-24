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
module spamm_bounding_box_2d

  implicit none

  !> A bounding box.
  type :: bounding_box_2d

     !> Lower row bound.
     integer :: lower(0:1) = -1

     !> Upper row bound.
     integer :: upper(0:1) = -1

   contains

     procedure :: to_string => bounding_box_2d_to_string

  end type bounding_box_2d

contains

  !> String representation of bounding box.
  !!
  !! @param self The bounding box.
  !!
  !! @return The string representation.
  function bounding_box_2d_to_string (self) result(string)

    character(len = 1000) :: string
    class(bounding_box_2d), intent(in) :: self

    character(len = 200) :: temp

    write(temp, *) self%lower(0)
    write(string, "(A)") "bbox = [ "//trim(adjustl(temp))

    write(temp, *) self%upper(0)
    write(string, "(A)") trim(string)//":"//trim(adjustl(temp))

    write(temp, *) self%lower(1)
    write(string, "(A)") trim(string)//", "//trim(adjustl(temp))

    write(temp, *) self%upper(1)
    write(string, "(A)") trim(string)//":"//trim(adjustl(temp))//" ]"

  end function bounding_box_2d_to_string

end module spamm_bounding_box_2d
