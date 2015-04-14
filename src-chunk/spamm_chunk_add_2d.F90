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

!> Module for adding 2D chunks.
module spamm_chunk_add_2d

  use spamm_chunk_type_2d

  implicit none

contains

  !> Add two 2-d chunks.
  !!
  !! \f[ C \leftarrow A + B \f]
  !!
  !! @param A Chunk A.
  !! @param B Chunk B.
  !! @return The result.
  function chunk_add_2d_2d (A, B) result (C)

    use spamm_chunk_copy_2d

    type(chunk_2d_t), pointer :: C
    type(chunk_2d_t), pointer, intent(in) :: A, B

    if(associated(C)) then
       deallocate(C)
    endif

    if(.not. associated(A) .and. .not. associated(B)) then
       return
    elseif(associated(A) .and. .not. associated(B)) then
       call chunk_copy_2d(A, C)
       return
    endif

  end function chunk_add_2d_2d

end module spamm_chunk_add_2d
