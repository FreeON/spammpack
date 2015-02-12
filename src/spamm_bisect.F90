!> Bisection functions.
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
module spamm_bisect

#include "spamm_utility_macros.h"

  implicit none

contains

  !> Bisect an index range.
  !!
  !! The index range is inclusive, i.e. all integers in the interval
  !! \f$ [ i_{l}, i_{u} ] \f$ are included. The function splits this
  !! interval into two. The split intervals are \f$ [ i_{l}, i_{m} ]
  !! \f$ and \f$ [ i_{m}+1, i_{u} ] \f$.
  !!
  !! @param lower The lower index.
  !! @param upper The upper index.
  !!
  !! @return The middle.
  function bisect (lower, upper) result (middle)

    use spamm_globals
    use spamm_strings

    integer :: middle
    integer, intent(in) :: lower, upper

    integer :: width
    character(len=1000) :: error_message

    width = upper-lower+1
    if(width > SPAMM_BLOCK_SIZE) then
       middle = width/2+lower-1
    else if(width == SPAMM_BLOCK_SIZE) then
       middle = -1
    else
       error_message = "can not bisect ["//trim(to_string(lower)) &
            //", "//trim(to_string(upper))//"] with block size " &
            //trim(to_string(SPAMM_BLOCK_SIZE))
       LOG_FATAL(error_message)
       error stop
    endif

  end function bisect

end module spamm_bisect
