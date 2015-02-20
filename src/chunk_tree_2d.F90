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
module chunk_tree_2d

#include "spamm_utility_macros.h"

  use spamm_decoration_2d
  use spamm_extra_2d
  use spamm_globals

  implicit none

  !> The tree-node type.
  type :: chunk_node_2d

     !> The tree node decoration.
     type(decoration_2d) :: decoration

     !> Some extra information.
     type(extra_2d) :: extra

  end type chunk_node_2d

  !> The leaf-node matrix type.
  type :: chunk_data_2d

     !> The matrix.
     real(SPAMM_KIND) :: data(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE)

  end type chunk_data_2d

  !> The 2D chunk type.
  type :: chunk_2d

     !> The tree nodes.
     type(chunk_node_2d), allocatable :: node(:)

     !> The matrices at the leaves.
     type(chunk_data_2d), allocatable :: data(:, :)

   contains

     procedure :: to_string => chunk_2d_to_string

  end type chunk_2d

contains

  function new_chunk_2d (N) result (chunk)

    type(chunk_2d), pointer :: chunk
    integer, intent(in) :: N
    integer :: depth, number_nodes

    allocate(chunk)

    depth = 0
    number_nodes = 0
    do while(.true.)
       if(SPAMM_BLOCK_SIZE*2**depth >= N) then
          exit
       endif
       number_nodes = number_nodes+4**depth
       depth = depth+1
    enddo

    allocate(chunk%node(number_nodes))
    allocate(chunk%data(2**depth, 2**depth))

  end function new_chunk_2d

  function chunk_2d_to_string (self) result (string)

    use spamm_strings

    character(len=1000) :: string
    class(chunk_2d), intent(in) :: self

    string = "chunk:"
    write(string, "(A)") trim(string)//" "//trim(to_string(size(self%node)))//" nodes"
    write(string, "(A)") trim(string)//", self: "//trim(to_string(storage_size(self)/8))//" bytes"
    write(string, "(A)") trim(string)//", self%node: "//trim(to_string(storage_size(self%node)/8))//" bytes"
    write(string, "(A)") trim(string)//", self%data: "//trim(to_string(storage_size(self%data)/8))//" bytes"

  end function chunk_2d_to_string

end module chunk_tree_2d
