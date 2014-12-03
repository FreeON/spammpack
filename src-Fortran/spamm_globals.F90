!> Defines global parameters and variables.
!!
!! @todo Move matrix related global parameters or values into appropriate types
!! and make things much less global.
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
!! @author Nicolas Bock nicolasbock@freeon.org
MODULE spamm_globals

  use spamm_types

  IMPLICIT NONE

  !> The size of the basic submatrix blocks.
  integer, parameter :: SPAMM_BLOCK_SIZE = 4

  !> The norm cutoff for tasked recursion.
  real(spamm_kind), parameter :: spamm_recursion_normd_cutoff = 1e-4

  !> @deprecated The depth of the matrix tree.
  INTEGER :: SpAMM_TOTAL_DEPTH

  !> @deprecated The size of the unpadded matrix.
  INTEGER :: SpAMM_MATRIX_DIMENSION

  !> @deprecated The size of the padded matrix.
  INTEGER :: SpAMM_PADDED_MATRIX_DIMENSION

  !> Cutoff the tree depth at some predefined maximum depth.
  INTEGER :: SpAMM_RECURSION_DEPTH_CUTOFF

END MODULE spamm_globals
