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

!> Module containing the 2D chunk data types.
!!
!! \todo
!!    - Implement setting a matrix from dense.
!!        - Update norms.
!!        - Check/Update bounds.
module spamm_chunk_type_2d

  implicit none

  !> The tree-node type.
  type :: chunk_node_2d_t
     sequence
     !> The number of non-zero elements.
     double precision :: number_nonzeros = 0
     !> The square of the Frobenius norm.
     double precision :: norm2 = 0
  end type chunk_node_2d_t

  !> The leaf-node matrix type.
  type :: chunk_data_2d_t
     sequence
     !> The matrix.
     double precision :: data(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE) = 0
  end type chunk_data_2d_t

  !> Bounds.
  type :: chunk_bounds_t
     sequence
     !> The lower and upper bound along a dimension.
     integer :: bounds(0:1) = [ -1, -1 ]
  end type chunk_bounds_t

  !> The 2D chunk type.
  !!
  !! The components are ordered such that the matrix data comes
  !! first. The `sequence`but attribute guarantees that the storage order
  !! of the components is the same as layed out in the type
  !! definition. All components are guaranteed to be allocated
  !! contiguously (Section 5.3.7 F08 Standard).
  type :: chunk_2d_t
     sequence
     !> The matrices at the leaves.
     type(chunk_data_2d_t) :: data(SPAMM_CHUNK_BLOCKS, SPAMM_CHUNK_BLOCKS)
     !> The tree nodes. The tree are stored in a complete quadtree,
     !! i.e.
     !!
     !! The nodes can be indexed from a starting index \f$ i \f$ as
     !! follows (this assumes that the root is at index 1):
     !!   - Children of node \f$ i \f$: \f$ \rightarrow 4 (i-1) + \{
     !!     2, 3, 4, 5 \} \f$.
     !!   - Parent of node \f$ i \f$: \f$ \lfloor \frac{ i-2 }{4}
     !!     \rfloor + 1 \f$.
     type(chunk_node_2d_t) :: node(SPAMM_CHUNK_NODES)
     !> Lower index bound.
     type(chunk_bounds_t) :: lower = chunk_bounds_t([ 1, 1 ])
     !> Upper index bound.
     type(chunk_bounds_t) :: upper = chunk_bounds_t([ SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE ])
  end type chunk_2d_t
  !> Interface to the equals functions.
  interface equals
     module procedure chunk_bounds_equals
  end interface equals

contains

  !> Test for equality of chunk_bounds_t
  !!
  !! @param A Bounds A.
  !! @param B Bounds B.
  !! @return .True. if the bounds are equal to each other.
  function chunk_bounds_equals (A, B) result(are_equal)

    type(chunk_bounds_t), intent(in) :: A
    type(chunk_bounds_t), intent(in) :: B
    logical :: are_equal

    if(A%bounds(0) == B%bounds(0) .and. A%bounds(1) == B%bounds(1)) then
       are_equal = .true.
    else
       are_equal = .false.
    endif

  end function chunk_bounds_equals

  !> Convert the meta-data of a chunk to a string.
  !!
  !! @param A The chunk.
  !! @return The string representation.
  function chunk_2d_to_string (A) result (string)

    use spamm_chunk_strings

    character(len=1000) :: string
    type(chunk_2d_t), intent(in) :: A

    string = "chunk:"
    write(string, "(A)") trim(string)//" N_c: "//trim(to_string(SPAMM_CHUNK_SIZE))
    write(string, "(A)") trim(string)//", box: [["//trim(to_string(A%lower%bounds(0)))// &
         ", "//trim(to_string(A%upper%bounds(0)))//"]"
    write(string, "(A)") trim(string)//", ["//trim(to_string(A%lower%bounds(1)))// &
         ", "//trim(to_string(A%upper%bounds(1)))//"]]"
    write(string, "(A)") trim(string)//"; N_b: "//trim(to_string(SPAMM_BLOCK_SIZE))
    write(string, "(A)") trim(string)//"; "//trim(to_string(size(A%node)))//" nodes"
    write(string, "(A)") trim(string)//"; "//trim(to_string(size(A%data, 1)**2))//" leaves"
    write(string, "(A)") trim(string)//"; total: "//trim(to_string(storage_size(A)/8))//" B"
    write(string, "(A)") trim(string)//"; node: "//trim(to_string(storage_size(A%node)/8))//" B"
    write(string, "(A)") trim(string)//"; basic: "// &
         trim(to_string(storage_size(A%data)/8))//" B"

  end function chunk_2d_to_string

  !> Write memory layout of chunk to string.
  !!
  !! @param A The chunk.
  !! @return The string representation.
  function chunk_2d_memory_layout (A) result (string)

    use spamm_chunk_strings
    use, intrinsic :: iso_C_binding

    character(len=1000) :: string
    type(chunk_2d_t), pointer, intent(in) :: A

    character(len=100) :: temp
    integer(c_intptr_t) :: ptr

    ptr = transfer(c_loc(A), ptr)+storage_size(A)

    write(temp, "(Z32)") ptr
    write(string, "(A)") "chunk layout:"//C_NEW_LINE// &
         trim(to_string(c_loc(A)))//": chunk start"//C_NEW_LINE// &
         trim(to_string(c_loc(A%data(1, 1))))//": data(1, 1)"//C_NEW_LINE// &
         trim(to_string(c_loc(A%data(2, 1))))//": data(2, 1)"//C_NEW_LINE// &
         trim(to_string(c_loc(A%data(1, 2))))//": data(1, 2)"//C_NEW_LINE// &
         trim(to_string(c_loc(A%data(SPAMM_CHUNK_BLOCKS, 1))))// &
         ": data("//trim(to_string(SPAMM_CHUNK_BLOCKS))//", 1)"//C_NEW_LINE// &
         trim(to_string(c_loc(A%data(SPAMM_CHUNK_BLOCKS, SPAMM_CHUNK_BLOCKS))))// &
         ": data("//trim(to_string(SPAMM_CHUNK_BLOCKS))//", "// &
         trim(to_string(SPAMM_CHUNK_BLOCKS))//")"//C_NEW_LINE// &
         trim(to_string(c_loc(A%node(1))))//": node(1)"//C_NEW_LINE// &
         trim(to_string(c_loc(A%node(SPAMM_CHUNK_NODES))))// &
         ": node("//trim(to_string(SPAMM_CHUNK_NODES))//")"//C_NEW_LINE// &
         trim(to_string(c_loc(A%lower%bounds(1))))//": lower(1)"//C_NEW_LINE// &
         trim(to_string(c_loc(A%upper%bounds(1))))//": upper(1)"//C_NEW_LINE// &
         "0x"//trim(adjustl(temp))//": chunk end"

  end function chunk_2d_memory_layout

  !> Convert a chunk to a dense matrix.
  !!
  !! @param A The chunk.
  !! @return The dense matrix.
  function chunk_2d_to_dense (A) result (B)

    double precision :: B(SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE)
    type(chunk_2d_t), pointer, intent(in) :: A

    integer :: i, j, k, l

    do i = 1, SPAMM_CHUNK_BLOCKS
       do j = 1, SPAMM_CHUNK_BLOCKS
          k = (i-1)*SPAMM_BLOCK_SIZE+1
          l = (j-1)*SPAMM_BLOCK_SIZE+1
          B(k:k+SPAMM_BLOCK_SIZE-1, l:l+SPAMM_BLOCK_SIZE-1) = A%data(i, j)%data
       enddo
    enddo

  end function chunk_2d_to_dense

  !> Get a matrix element from a chunk. The matrix indices are based
  !! off of the stored bounding box, i.e. may not start at 1.
  !!
  !! @param i The row index.
  !! @param j The column index.
  !! @param A The chunk.
  !! @return The matrix element.
  function chunk_2d_get (i, j, A) result(Aij)

    double precision :: Aij
    integer, intent(in) :: i, j
    type(chunk_2d_t), intent(in) :: A

    integer :: i_leaf(0:1), i_block(0:1)

    Aij = 0

    if(i < A%lower%bounds(0) .or. i > A%upper%bounds(0) .or. &
         j < A%lower%bounds(1) .or. j > A%upper%bounds(1)) return

    i_leaf(0) = (i-A%lower%bounds(0))/SPAMM_BLOCK_SIZE+1
    i_leaf(1) = (j-A%lower%bounds(1))/SPAMM_BLOCK_SIZE+1

    i_block(0) = mod(i-A%lower%bounds(0), SPAMM_BLOCK_SIZE)+1
    i_block(1) = mod(j-A%lower%bounds(1), SPAMM_BLOCK_SIZE)+1

    Aij = A%data(i_leaf(0), i_leaf(1))%data(i_block(0), i_block(1))

  end function chunk_2d_get

  !> Calculate the matrix norm of a dense matrix. The norm is the
  !! Frobenius norm.
  !!
  !! @param A The dense matrix.
  !! @return The matrix norm.
  function matrix_norm (A) result(norm)

    double precision, intent(in) :: A(:, :)
    double precision :: norm

    norm = 0

  end function matrix_norm

  !> Check the internal consistency and correctness of a chunk. For
  !! instance, the function verifies that all norms are correct.
  !!
  !! In detail:
  !!   - Verify all norms.
  !!   - Verify non-zero counts.
  !!
  !! @param A The chunk.
  !! @return Whether the chunk passed the check or not.
  function chunk_2d_check (A) result(is_verified)

    type(chunk_2d_t), pointer, intent(in) :: A
    logical :: is_verified

    integer :: i, j
    double precision :: norms(SPAMM_CHUNK_BLOCKS, SPAMM_CHUNK_BLOCKS)

    is_verified = .true.

    if(.not. associated(A)) return
    if(A%lower%bounds(0) <= A%upper%bounds(0)) return
    if(A%lower%bounds(1) <= A%upper%bounds(1)) return

    do i = 1, SPAMM_CHUNK_BLOCKS
       do j = 1, SPAMM_CHUNK_BLOCKS
          norms(i, j) = matrix_norm(A%data(i, j)%data)
       enddo
    enddo

  end function chunk_2d_check

end module spamm_chunk_type_2d
