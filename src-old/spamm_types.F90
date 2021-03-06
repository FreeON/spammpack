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
module spamm_types

#ifdef _OPENMP
  use omp_lib
#endif

  use spamm_real_precision

  implicit none

  !> Define the number zero.
  !real(SPAMM_KIND), parameter :: SpAMM_Zero = 0D0

  !> Define the number 1/2.
  !real(SPAMM_KIND), parameter :: SpAMM_Half = 5D-1

  !> Define the number 1.
  real(SPAMM_KIND), parameter :: SpAMM_One = 1D0

  !> Define the number 2.
  real(SPAMM_KIND), parameter :: SpAMM_Two = 2D0

  !> Define the number 4.
  real(SPAMM_KIND), parameter :: SpAMM_Four = 4D0

  !> Define the number 8.
  real(SPAMM_KIND), parameter :: SpAMM_Eight = 8D0

  !> Bigest machine double for ONX_KIND
  !real(SPAMM_KIND), parameter :: SpAMM_BIG_DBL = huge(SpAMM_One)

  !> Bigest machine int for int*4
  integer, parameter :: SpAMM_BIG_INT = 2**28

  type :: decoration_1d
     integer                     :: N = -1
     integer                     :: N_padded = -1
     integer                     :: depth = -1
     integer, dimension(1:2)     :: bb = -1
     real(kind(0d0))             :: number_nonzeros = -1
     real(SPAMM_KIND)            :: norm2 = -1
  end type decoration_1d

  type :: decoration_2d
     integer                     :: N = -1
     integer                     :: N_padded = -1
     integer                     :: depth = -1
     integer, dimension(1:2,0:1) :: bb = -1
     real(kind(0d0))             :: number_nonzeros = -1
     real(SPAMM_KIND)            :: norm2 = -1
  end type decoration_2d

  type :: SpAMM_tree_1d
     type(decoration_1d)    :: decoration 
     type(SpAMM_tree_1d), pointer :: child_0
     type(SpAMM_tree_1d), pointer :: child_1
     real(SPAMM_KIND), dimension(:), allocatable :: vect
  end type SpAMM_tree_1d  

  type :: SpAMM_tree_2d
     type(decoration_2d)    :: decoration
     type(SpAMM_tree_2d), pointer :: child_00 => null()
     type(SpAMM_tree_2d), pointer :: child_01 => null()
     type(SpAMM_tree_2d), pointer :: child_11 => null()
     real(SPAMM_KIND), dimension(:), allocatable :: blok
  end type SpAMM_tree_2d

  type :: SpAMM_tree_2d_full
     type(decoration_2d)    :: decoration
     type(SpAMM_tree_2d), pointer :: child_00 => null()
     type(SpAMM_tree_2d), pointer :: child_01 => null()
     type(SpAMM_tree_2d), pointer :: child_10 => null()
     type(SpAMM_tree_2d), pointer :: child_11 => null()
     real(SPAMM_KIND), dimension(:), allocatable :: blok
  end type SpAMM_tree_2d_full

  !> A vector type.
  type spamm_matrix_order_1

     !> The number of entries.
     integer :: N = -1

     !> The padded vector dimension.
     integer :: N_padded

     !> The tree depth. The root tier is 0, tier == depth is the leaf node
!! tier, i.e. the tier at which actual matrix elements are stored.
     integer :: depth = -1

     !> The Frobenius norm.
     real(SPAMM_KIND) :: norm = 0

     !> The root of the binary tree.
     type(bitree), pointer :: root => null()

     !> The number of non-zero elements.
     real(kind(0d0)) :: number_nonzeros = 0

  end type spamm_matrix_order_1

  !> Binary tree data structure.
  TYPE BiTree

     !> The norm.
     REAL(SPAMM_KIND) :: Norm = 0

     !> The number of non-zero elements.
     real(kind(0d0)) :: number_nonzeros = 0

     !> The lower row index.
     integer :: i_lower = -1

     !> The upper row index.
     integer :: i_upper = -1

     !> The pointer to the left bisecting subtree.
     TYPE(BiTree), POINTER :: Sect1 => NULL()

     !> The pointer to the right bisecting subtree.
     TYPE(BiTree), POINTER :: Sect2 => NULL()

     !> The vector data.
     REAL(SPAMM_KIND), DIMENSION(:), ALLOCATABLE :: Vect

  END TYPE BiTree

  !> Matrix (2nd order) type.
  type spamm_matrix_order_2

     !> The number or rows.
     integer :: M = -1

     !> The number of columns.
     integer :: N = -1

     !> The tree depth. The root tier is 0, tier == depth is the leaf node
!! tier, i.e. the tier at which actual matrix elements are stored.
     integer :: depth = -1

     !> The Frobenius norm.
     REAL(SPAMM_KIND) :: norm = 0

     !> The padded matrix dimension. Since the padded matrix is always square,
!! we only store one number here.
     integer :: N_padded

     !> The root quadtree pointer.
     type(qutree), pointer :: root => null()

     !> The number of non-zero elements.
     REAL(kind(0d0)) :: number_nonzeros = 0

     !> The number of operations (updated by a Multiply), i.e. the number of
!! dense matrix products of size spamm_globals::spamm_block_size x
!! spamm_globals::spamm_block_size.
     REAL(kind(0d0)) :: number_operations = 0

  end type spamm_matrix_order_2

  !> Quaternary tree data structure.
  type QuTree

     !> The lower row index.
     integer :: i_lower = -1

     !> The upper row index.
     integer :: i_upper = -1

     !> The lower column index.
     integer :: j_lower = -1

     !> The upper column index.
     integer :: j_upper = -1

     !> The Frobenious norm.
     REAL(SPAMM_KIND) :: Norm = 0

     !> The pointer to the subtree in quadrant 11.
     TYPE(QuTree), POINTER :: Quad11 => NULL()

     !> The pointer to the subtree in quadrant 12.
     TYPE(QuTree), POINTER :: Quad12 => NULL()

     !> The pointer to the subtree in quadrant 21.
     TYPE(QuTree), POINTER :: Quad21 => NULL()

     !> The pointer to the subtree in quadrant 22.
     TYPE(QuTree), POINTER :: Quad22 => NULL()

     !> The matrix data.
     REAL(SPAMM_KIND), DIMENSION(:, :), ALLOCATABLE :: Blok

#ifdef SPAMM_STORE_TRANSPOSE
     !> The transposed block.
     real(SPAMM_KIND), dimension(:, :), ALLOCATABLE :: transpose_block
#endif

     !> The number of non-zero elements.
     REAL(Kind(0d0)) :: number_nonzeros = 0

     !> The number of operations (updated by a Multiply), i.e. the number of
!! dense matrix products of size spamm_globals::spamm_block_size x
!! spamm_globals::spamm_block_size.
     REAL(Kind(0d0)) :: number_operations = 0

#ifdef _OPENMP
     !> Block OpenMP lock
     integer(kind = OMP_LOCK_KIND) :: lock
#endif

  end type QuTree

end module spamm_types
