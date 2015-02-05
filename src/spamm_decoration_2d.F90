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

  use spamm_real_precision

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

     !> The square of the Frobenius norm.
     real(SPAMM_KIND) :: norm2 = -1

     !> The number of non-zero elements.
     real(kind(0d0)) :: number_nonzeros = -1

   contains

     procedure :: to_string => decoration_2d_to_string

  end type decoration_2d

contains

  !> String representation of decoration.
  !!
  !! @param self The node decoration.
  !!
  !! @return The string representation.
  function decoration_2d_to_string (self) result (string)

    character(len = 1000) :: string
    class(decoration_2d), intent(in) :: self

    character(len = 200) :: temp

    write(temp, *) self%N
    write(string, "(A)") "N = "//trim(adjustl(temp))

    write(temp, *) self%N_padded
    write(string, "(A)") trim(string)//", N_padded = "//trim(adjustl(temp))

    write(temp, *) self%depth
    write(string, "(A)") trim(string)//", depth = "//trim(adjustl(temp))

    write(temp, "(ES15.5)") self%norm2
    write(string, "(A)") trim(string)//", norm2 = "//trim(adjustl(temp))

    write(temp, "(ES15.5)") self%number_nonzeros
    write(string, "(A)") trim(string)//", nnonz = "//trim(adjustl(temp))

  end function decoration_2d_to_string

  !> Merge two decorations.
  !!
  !! @param A Decoration A.
  !! @param B Decoration B.
  function decoration_2d_merge (A, B) result (C)

    type(decoration_2d) :: C
    type(decoration_2d), intent(in) :: A, B

    if(A%N /= B%N) then
       LOG_FATAL("dimension mismatch in N")
       error stop
    endif

    if(A%N_padded /= B%N_padded) then
       LOG_FATAL("padded dimension mismatch")
       error stop
    endif

    if(A%depth /= B%depth) then
       LOG_FATAL("depth mismatch")
       error stop
    endif

    C%N = A%N
    C%N_padded = A%N_padded
    C%depth = A%depth
    C%norm2 = A%norm2+B%norm2
    C%number_nonzeros = A%number_nonzeros+B%number_nonzeros

  end function decoration_2d_merge

end module spamm_decoration_2d
