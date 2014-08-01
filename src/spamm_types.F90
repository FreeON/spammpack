!> Defines derived types used in SpAMMPACK.
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
module spamm_types

  !$ USE OMP_LIB

  IMPLICIT NONE

  !> Define integer of length 2.
  INTEGER, PARAMETER :: INT1 = SELECTED_INT_KIND(2)  !--Integer*1

  !> Define integer of length 2.
  INTEGER, PARAMETER :: INT2 = SELECTED_INT_KIND(4)  !--Integer*2

  !> Define integer of length 4.
  INTEGER, PARAMETER :: INT4 = SELECTED_INT_KIND(9)  !--Integer*4

  !> Define integer of length 8.
  INTEGER, PARAMETER :: INT8 = SELECTED_INT_KIND(18) !--Integer*8

  !> Define float of length 4.
  INTEGER, PARAMETER :: SpAMM_SINGLE = KIND(0.0e0)   !--Real*4

  !> Define float of length 8.
  INTEGER, PARAMETER :: SpAMM_DOUBLE = KIND(0.0d0)   !--Real*8

#ifdef SPAMM_SINGLE
  !> Define floating point type to single.
  INTEGER, PARAMETER :: SpAMM_KIND = spamm_single
#else
  !> Define floating point type to double.
  INTEGER, PARAMETER :: SpAMM_KIND = spamm_double
#endif

  !> Define the number zero.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Zero = 0D0

  !> Define the number 1/2.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Half = 5D-1

  !> Define the number 1.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_One = 1D0

  !> Define the number 2.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Two = 2D0

  !> Define the number 4.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Four = 4D0

  !> Define the number 8.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Eight = 8D0

  !> Bigest machine double for ONX_KIND
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_BIG_DBL = HUGE(SpAMM_One)

  !> Bigest machine int for int*4
  INTEGER, PARAMETER :: SpAMM_BIG_INT = 2**28

  !> Binary tree data structure.
  TYPE BiTree
    !> The norm.
    REAL(SpAMM_KIND) :: Norm

    !> The pointer to the left bisecting subtree.
    TYPE(BiTree), POINTER :: Sect0 => NULL()

    !> The pointer to the right bisecting subtree.
    TYPE(BiTree), POINTER :: Sect1 => NULL()

    !> The vector data.
    REAL(SpAMM_KIND), DIMENSION(:), ALLOCATABLE :: Vect
  END TYPE BiTree

  !> Matrix (2nd order) type.
  type spamm_matrix_2nd_order
    !> The root quadtree pointer.
    type(qutree), pointer :: root => null()

    !> The number of non-zero elements.
    REAL(SpAMM_DOUBLE) :: number_nonzeros

    !> The number of operations (updated by a Multiply), i.e. the number of dense matrix products of size
    !! spamm_globals::spamm_block_size x spamm_globals::spamm_block_size.
    REAL(SpAMM_DOUBLE) :: number_operations

  end type spamm_matrix_2nd_order

  !> Quaternary tree data structure.
  TYPE QuTree
    !> The Frobenious norm.
    REAL(SpAMM_KIND) :: Norm = 0

    !> The pointer to the subtree in quadrant 11.
    TYPE(QuTree), POINTER :: Quad11 => NULL()

    !> The pointer to the subtree in quadrant 12.
    TYPE(QuTree), POINTER :: Quad12 => NULL()

    !> The pointer to the subtree in quadrant 21.
    TYPE(QuTree), POINTER :: Quad21 => NULL()

    !> The pointer to the subtree in quadrant 22.
    TYPE(QuTree), POINTER :: Quad22 => NULL()

    !> The matrix data.
    REAL(SpAMM_KIND), DIMENSION(:, :), ALLOCATABLE :: Blok

    !> The number of non-zero elements.
    REAL(SpAMM_DOUBLE) :: number_nonzeros

    !> The number of operations (updated by a Multiply), i.e. the number of dense matrix products of size
    !! spamm_globals::spamm_block_size x spamm_globals::spamm_block_size.
    REAL(SpAMM_DOUBLE) :: number_operations

#ifdef _OPENMP
    !> Block lock
    INTEGER(KIND = OMP_LOCK_KIND) :: lock
#endif
  END TYPE QuTree

  !> A type for performance measurements.
  TYPE Stats
    !> The time.
    REAL(SpAMM_DOUBLE) :: Time

    !> Some count.
    INTEGER            :: Count

    !> The name of a function.
    CHARACTER(LEN=50)  :: Routine
  END TYPE Stats

CONTAINS

end module spamm_types
