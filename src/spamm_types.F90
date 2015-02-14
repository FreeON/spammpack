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
!
!
!----------------------------------------------------------------------------------
! The "long weekend hack".  Going for a functional implementation
! that fits tightly in the edit window (mc,2/2015).

MODULE spamm_types

  USE  spamm_real_precision
  USE  spamm_numbers

  ! SpAMM Garnishments ________________ SGARN _________________

  ! Enrichments of the tree_1d ...
  type :: decoration_1d
     !> Integer dimension of the native (non-padded) vector
     integer                               :: NDimn
     !> Axis-aligned brakets for the (i-leftmost) index space
     integer,  pointer, dimension(0:1)     :: BndBx
     !> Square of the F-norm.
     real(SPAMM_KIND)                      :: Norm2 = -1
     !> Float Ops needed accumulated to this level 
     real(kind(0d0))                       :: FlOps = -1
     !> The number of non-zero elements to this level
     real(kind(0d0))                       :: Non0s = -1
  end type decoration_1d

  ! Enrichments of the tree_2d ...
  type :: decoration_2d
     !> Integer dimension of the native (non-padded) matrix 
     integer,           dimension(1:2)     :: NDimn
     !> Axis-aligned bounding box for the (i-j) index space
     integer,  pointer, dimension(0:1,1:2) :: BndBx
     !> Square of the F-norm.
     real(SPAMM_KIND)                      :: Norm2 = -1
     !> Float Ops needed accumulated to this level 
     real(kind(0d0))                       :: FlOps = -1
     !> The number of non-zero elements to this level 
     real(kind(0d0))                       :: Non0s = -1
  end type decoration_2d

  ! SpAMM algebraic data types _______________ SADTS ____________________

  ! The tree_1d (vector) type ...
  type :: SpAMM_tree_1d
     type(SpAMM_decoration_1d)             :: frill
     type(SpAMM_tree_1d_symm), pointer     :: child_0 => null()
     type(SpAMM_tree_1d_symm), pointer     :: child_1 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_1d

  ! The tree_2d matrix types...

  ! full:
  type :: SpAMM_tree_2d_full
     type(SpAMM_decoration_2d)             :: frill
     type(SpAMM_tree_2d_full), pointer     :: child_00 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_01 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_10 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_11 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_2d_full

  ! symmetric (SPD/Hermetian) :
  type :: SpAMM_tree_2d_symm
     type(SpAMM_decoration_2d)             :: frill
     type(SpAMM_tree_2d_symm), pointer     :: child_00 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_01 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_11 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_2d_symm

  ! --
  implicit none

  ! --
contains

  ! --
end module spamm_types ! ... and we're out ...
