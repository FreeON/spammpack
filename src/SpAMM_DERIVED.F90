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
MODULE SpAMM_DERIVED
  !$ USE OMP_LIB
  IMPLICIT NONE
  PUBLIC
  !> Define integer of length 2.
  INTEGER, PARAMETER :: INT1 = SELECTED_INT_KIND(2)  !--Integer*1

  !> Define integer of length 2.
  INTEGER, PARAMETER :: INT2 = SELECTED_INT_KIND(4)  !--Integer*2
  !> Define integer of length 4.
  INTEGER, PARAMETER :: INT4 = SELECTED_INT_KIND(9)  !--Integer*4
  !> Define integer of length 8.
  INTEGER, PARAMETER :: INT8 = SELECTED_INT_KIND(18) !--Integer*8
  !> Define float of length 4.
  INTEGER, PARAMETER :: SpAMM_SINGLE = KIND(0.0E0)   !--Real*4
  !> Define float of length 8.
  INTEGER, PARAMETER :: SpAMM_DOUBLE = KIND(0.0D0)   !--Real*8
  !> Define the float type used for SpAMM.
#ifdef SPAMM_SINGLE
  INTEGER, PARAMETER :: SpAMM_KIND = SpAMM_SINGLE
#else
  INTEGER, PARAMETER :: SpAMM_KIND = SpAMM_DOUBLE
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
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_BIG_DBL=HUGE(SpAMM_One)           
  !> Bigest machine int for int*4
  INTEGER,      PARAMETER     :: SpAMM_BIG_INT=2**28  
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
  ! --------------------------------------------------
  ! Quad structures
  ! --------------------------------------------------
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
#ifdef _OPENMP
    !> Block lock 
    INTEGER(KIND = OMP_LOCK_KIND) :: Lock
#endif
  END TYPE QuTree

  TYPE Stats
     REAL(SpAMM_DOUBLE) :: Time
     INTEGER            :: Count
     CHARACTER(LEN=50)  :: Routine
  END TYPE Stats

  !! This type includes several different norms.
  TYPE SpAMM_Norm

    !> The Frobenius norm.
    !!
    !! @f$ \Vert A \Vert_{F} = @f$.
    REAL(SpAMM_KIND) :: FrobeniusNorm = SpAMM_Zero

    !> The max norm.
    !!
    !! @f$ \max_{ij} A_{ij} @f$.
    REAL(SpAMM_KIND) :: MaxNorm = SpAMM_Zero

  END TYPE SpAMM_Norm

CONTAINS

END MODULE SpAMM_DERIVED

!!$
!!$
!!$  FUNCTION QuTreeConstructor () RESULT(q)
!!$
!!$    TYPE(QuTree) :: q
!!$
!!$    q%Quad11 => NULL()
!!$    q%Quad12 => NULL()
!!$    q%Quad21 => NULL()
!!$    q%Quad22 => NULL()
!!$
!!$  END FUNCTION QuTreeConstructor
!!$
!!$  TYPE QuLeaf
!!$
!!$    REAL(SpAMM_KIND), DIMENSION(16, 16) :: Blok
!!$
!!$  END TYPE QuLeaf
!!$
!!$  !INTERFACE QuTree
!!$  !  MODULE PROCEDURE QuTreeConstructor
!!$  !END INTERFACE QuTree
!!$
!!$  TYPE QuLink
!!$     TYPE(QuLink), POINTER  :: Next
!!$     TYPE(QuTree), POINTER  :: Quad
!!$     REAL(SpAMM_KIND)       :: Wght
!!$     INTEGER                :: Hash
!!$     INTEGER                :: i
!!$     INTEGER, DIMENSION(2,2):: Box
!!$  END TYPE QuLink
!!$
!!$  ! QuTree container
!!$  TYPE QuTptr
!!$     TYPE(QuTree),POINTER :: T
!!$  END TYPE QuTptr
!!$
!!$  ! QuLink container
!!$  TYPE QuLptr
!!$     TYPE(QuLink),POINTER :: L
!!$  END TYPE QuLptr
!!$
!!$  ! --------------------------------------------------
!!$  ! Cube structures
!!$  ! --------------------------------------------------
!!$  TYPE CuTree
!!$     TYPE(CuTree), POINTER  :: Next
!!$     TYPE(QuTree), POINTER  :: qA
!!$     TYPE(QuTree), POINTER  :: qB
!!$     TYPE(QuTree), POINTER  :: qC
!!$     INTEGER, DIMENSION(3,2):: Box
!!$     INTEGER                :: Hash
!!$  END TYPE CuTree
!!$
!!$  TYPE CuLink
!!$     TYPE(CuLink), POINTER  :: Next
!!$     TYPE(QuTree), POINTER  :: QuadA
!!$     TYPE(QuTree), POINTER  :: QuadB
!!$     INTEGER, DIMENSION(3,2):: Box
!!$     REAL(SpAMM_KIND)       :: Wght
!!$     INTEGER                :: Hash
!!$     INTEGER                :: i
!!$  END TYPE CuLink
!!$
!!$  ! CuTree container
!!$  TYPE CuTPtr
!!$     TYPE(CuTree),POINTER :: T
!!$  END TYPE CuTPtr
!!$
!!$  ! CuLink container
!!$  TYPE CuLPtr
!!$     TYPE(CuLink),POINTER :: L
!!$  END TYPE CuLPtr
!!$ NULL()
!!$
!!$  END FUNCTION QuTreeConstructor
!!$
!!$END MODULE SpAMM_DERIVED
!!$
!!$
!!$  TYPE QuLeaf
!!$
!!$    REAL(SpAMM_KIND), DIMENSION(16, 16) :: Blok
!!$
!!$  END TYPE QuLeaf
!!$
!!$  !INTERFACE QuTree
!!$  !  MODULE PROCEDURE QuTreeConstructor
!!$  !END INTERFACE QuTree
!!$
!!$  TYPE QuLink
!!$     TYPE(QuLink), POINTER  :: Next
!!$     TYPE(QuTree), POINTER  :: Quad
!!$     REAL(SpAMM_KIND)       :: Wght
!!$     INTEGER                :: Hash
!!$     INTEGER                :: i
!!$     INTEGER, DIMENSION(2,2):: Box
!!$  END TYPE QuLink
!!$
!!$  ! QuTree container
!!$  TYPE QuTptr
!!$     TYPE(QuTree),POINTER :: T
!!$  END TYPE QuTptr
!!$
!!$  ! QuLink container
!!$  TYPE QuLptr
!!$     TYPE(QuLink),POINTER :: L
!!$  END TYPE QuLptr
!!$
!!$  ! --------------------------------------------------
!!$  ! Cube structures
!!$  ! --------------------------------------------------
!!$  TYPE CuTree
!!$     TYPE(CuTree), POINTER  :: Next
!!$     TYPE(QuTree), POINTER  :: qA
!!$     TYPE(QuTree), POINTER  :: qB
!!$     TYPE(QuTree), POINTER  :: qC
!!$     INTEGER, DIMENSION(3,2):: Box
!!$     INTEGER                :: Hash
!!$  END TYPE CuTree
!!$
!!$  TYPE CuLink
!!$     TYPE(CuLink), POINTER  :: Next
!!$     TYPE(QuTree), POINTER  :: QuadA
!!$     TYPE(QuTree), POINTER  :: QuadB
!!$     INTEGER, DIMENSION(3,2):: Box
!!$     REAL(SpAMM_KIND)       :: Wght
!!$     INTEGER                :: Hash
!!$     INTEGER                :: i
!!$  END TYPE CuLink
!!$
!!$  ! CuTree container
!!$  TYPE CuTPtr
!!$     TYPE(CuTree),POINTER :: T
!!$  END TYPE CuTPtr
!!$
!!$  ! CuLink container
!!$  TYPE CuLPtr
!!$     TYPE(CuLink),POINTER :: L
!!$  END TYPE CuLPtr
!!$
!!$
