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
!! @author Nicolas Bock nicolasbock@freeon.org
module spamm_utilities

  implicit none

  integer, parameter :: FATAL = -1
  integer, parameter :: DEBUG = 2

  !> Interface to to_string functions.
  interface to_string
    module procedure int_to_string
    module procedure double_to_string
    module procedure bitree_to_string
    module procedure qutree_to_string
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

  !> Convert a bitree node to a string.
  !!
  !! @param q The bitree node.
  !!
  !! @return The string representation.
  function bitree_to_string (q) result(str_q)

    use spamm_types

    type(bitree), pointer, intent(in) :: q
    character(len = 200) :: temp
    character(len = :), allocatable :: str_q

    if(associated(q)) then
      write(temp, "(A)") "b["// &
        to_string(q%i_lower)//":"// &
        to_string(q%i_upper)//"]"

      if(allocated(q%vect)) then
        write(temp, "(A)") trim(temp)//", vect["// &
          int_to_string(size(q%vect, 1))//"]"
      else
        write(temp, "(A)") trim(temp)//", vect not allocated"
      endif

      if(associated(q%sect1)) then
        write(temp, "(A)") trim(temp)//", sect1"
      endif

      if(associated(q%sect2)) then
        write(temp, "(A)") trim(temp)//", sect2"
      endif

    else
      write(temp, "(A)") "q not associated"
    endif

    str_q = trim(adjustl(temp))

  end function bitree_to_string

  !> Convert a qutree node to a string.
  !!
  !! @param q The qutree node.
  !!
  !! @return The string representation.
  function qutree_to_string (q) result(str_q)

    use spamm_types

    type(qutree), pointer, intent(in) :: q
    character(len = 200) :: temp
    character(len = :), allocatable :: str_q

    if(associated(q)) then
      write(temp, "(A)") "q["// &
        to_string(q%i_lower)//":"// &
        to_string(q%i_upper)//","// &
        to_string(q%j_lower)//":"// &
        to_string(q%j_upper)//"]"

      if(allocated(q%blok)) then
        write(temp, "(A)") trim(temp)//", blok["// &
          int_to_string(size(q%blok, 1))//","// &
          int_to_string(size(q%blok, 2))//"]"
      else
        write(temp, "(A)") trim(temp)//", block not allocated"
      endif

      if(associated(q%quad11)) then
        write(temp, "(A)") trim(temp)//", quad11"
      endif

      if(associated(q%quad12)) then
        write(temp, "(A)") trim(temp)//", quad12"
      endif

      if(associated(q%quad21)) then
        write(temp, "(A)") trim(temp)//", quad21"
      endif

      if(associated(q%quad22)) then
        write(temp, "(A)") trim(temp)//", quad22"
      endif

    else
      write(temp, "(A)") "q not associated"
    endif

    str_q = trim(adjustl(temp))

  end function qutree_to_string

  !> Print a log message.
  !!
  !! @param level The message level. The global debug_level has to be greater or equal to the level.
  !! @param message The message to print.
  subroutine write_log (level, message)

    integer, intent(in) :: level
    character(len = *), intent(in) :: message

    if(SPAMM_DEBUG_LEVEL >= level) then
      if(level < 0) then
        write(*, "(A)") "[FATAL] "//message
        error stop
      else
        write(*, "(A)") message
      endif
    endif

  end subroutine write_log

end module spamm_utilities
