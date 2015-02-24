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
     real(kind(0d0)) :: data(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE)

  end type chunk_data_2d

  !> The 2D chunk type.
  type :: chunk_2d

     !> The tree nodes.
     type(chunk_node_2d) :: node(SPAMM_CHUNK_NODES)

     !> The matrices at the leaves.
     type(chunk_data_2d) :: data(SPAMM_BLOCKS, SPAMM_BLOCKS)

  end type chunk_2d

contains

  function chunk_2d_to_string (self) result (string)

    use spamm_strings

    character(len=1000) :: string
    type(chunk_2d), intent(in) :: self

    string = "chunk:"
    write(string, "(A)") trim(string)//" N_chunk: "//trim(to_string(SPAMM_CHUNK_SIZE))
    write(string, "(A)") trim(string)//", N_block: "//trim(to_string(SPAMM_BLOCK_SIZE))
    write(string, "(A)") trim(string)//", "//trim(to_string(size(self%node)))//" nodes"
    write(string, "(A)") trim(string)//", "//trim(to_string(size(self%data, 1)**2))//" leaves"
    write(string, "(A)") trim(string)//", total: "//trim(to_string(storage_size(self)/8))//" bytes"
    write(string, "(A)") trim(string)//", per node: "//trim(to_string(storage_size(self%node)/8))//" bytes"
    write(string, "(A)") trim(string)//", per matrix: "//trim(to_string(storage_size(self%data)/8))//" bytes"

  end function chunk_2d_to_string

  function chunk_2d_to_dense (A) result (B)

    real(kind(0d0)) :: B(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)
    type(chunk_2d), pointer, intent(in) :: A

    integer :: i, j, k, l

    do i = 1, SPAMM_BLOCKS
       do j = 1, SPAMM_BLOCKS
          k = (i-1)*SPAMM_BLOCK_SIZE+1
          l = (j-1)*SPAMM_BLOCK_SIZE+1
          B(k:k+SPAMM_BLOCK_SIZE-1, l:l+SPAMM_BLOCK_SIZE-1) = A%data(i, j)%data
       enddo
    enddo

  end function chunk_2d_to_dense

end module chunk_tree_2d
