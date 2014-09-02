!> Defines utilitiy functions.
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
!! @author Nicolas Bock nicolas.bock@freeon.org
module spamm_utilities

  implicit none

  integer, parameter :: FATAL = -1

  !> Interface to to_string functions.
  interface to_string
    module procedure int_to_string
    module procedure double_to_string
  end interface to_string

contains

  !> Convert an integer to a string.
  !!
  !! @param i The integer
  !!
  !! @return The string representation.
  function int_to_string (i) result(str_i)

    integer, intent(in) :: i

    character(len = 100) :: temp
    character(len = :), allocatable :: str_i

    write(temp, *) i
    str_i = trim(adjustl(temp))

  end function int_to_string

  !> Convert a double precision real to a string.
  !!
  !! @param x The real
  !!
  !! @return The string representation.
  function double_to_string (x) result(str_x)

    real(kind(0d0)), intent(in) :: x
    character(len = 100) :: temp
    character(len = :), allocatable :: str_x

    write(temp, "(ES16.8E3)") x
    str_x = trim(adjustl(temp))

  end function double_to_string

  !> Print a log message.
  !!
  !! @param level The message level. The global debug_level has to be greater or equal to the level.
  !! @param message An array of strings to print. Each string will be printed on a separate line.
  subroutine write_log (level, message)

    integer, intent(in) :: level
    character(len = *), intent(in) :: message

    if(DEBUG_LEVEL >= level) then
      if(level < 0) then
        write(*, "(A)") "[FATAL] "//message
        error stop
      else
        write(*, "(A)") message
      endif
    endif

  end subroutine write_log

end module spamm_utilities
