!> @copyright
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

!> Module containing string helper functions.
module spamm_chunk_strings

  implicit none

  !> Interface to to_string functions.
  interface to_string
    module procedure c_ptr_to_string
    module procedure double_array_to_string
    module procedure double_to_string
    module procedure int_to_string
    module procedure long_int_to_string
    module procedure single_to_string
  end interface to_string

contains

  !> Convert an integer to a string.
  !!
  !! @param i The integer.
  !! @return The string representation.
  function int_to_string (i) result (str_i)

    integer, intent(in) :: i

    character(len=100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len=:), allocatable :: str_i
#else
    character(len=100) :: str_i
#endif

    write(temp, *) i
    str_i = trim(adjustl(temp))

  end function int_to_string

  !> Convert an integer to a string.
  !!
  !! @param i The integer.
  !! @return The string representation.
  function long_int_to_string (i) result (str_i)

    integer(selected_int_kind(15)), intent(in) :: i

    character(len=100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len=:), allocatable :: str_i
#else
    character(len=100) :: str_i
#endif

    write(temp, *) i
    str_i = trim(adjustl(temp))

  end function long_int_to_string

  !> Convert a single precision real to a string.
  !!
  !! @param x The real
  !!
  !! @return The string representation.
  function single_to_string (x) result (str_x)

    real, intent(in) :: x
    character(len=100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len=:), allocatable :: str_x
#else
    character(len=100), allocatable :: str_x
#endif

    write(temp, "(ES16.8E3)") x
    str_x = trim(adjustl(temp))

  end function single_to_string

  !> Convert an array of double precision reals to a string.
  !!
  !! @param x The real array.
  !! @return The string representation.
  function double_array_to_string (x) result (str_x)

    double precision, intent(in) :: x(:)
    character(len=1000) :: temp
    character(len=100) :: format_string
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len=:), allocatable :: str_x
#else
    character(len=100), allocatable :: str_x
#endif
    integer :: i

    write(format_string, "(A,I4,A)") "(", size(x), "ES16.8E3)"
    write(temp, format_string) (x(i), i = 1, size(x))
    str_x = trim(adjustl(temp))

  end function double_array_to_string

  !> Convert a double precision real to a string.
  !!
  !! @param x The real.
  !! @return The string representation.
  function double_to_string (x) result (str_x)

    double precision, intent(in) :: x
    character(len=100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len=:), allocatable :: str_x
#else
    character(len=100), allocatable :: str_x
#endif

    write(temp, "(ES16.8E3)") x
    str_x = trim(adjustl(temp))

  end function double_to_string

  !> Convert a c_ptr to a string.
  !!
  !! @bug This is not portable apparently. It breaks under the Intel compiler.
  !!
  !! @param ptr The pointer.
  !! @return The string representation.
  function c_ptr_to_string (ptr) result (string)

    use, intrinsic :: iso_C_binding

    type(c_ptr), intent(in) :: ptr
    character(len=100) :: string

    !write(string, "(Z32)") ptr
    !write(string, "(A)") "0x"//trim(adjustl(string))

  end function c_ptr_to_string

end module spamm_chunk_strings
