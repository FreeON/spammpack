!> Defines operations on SpAMM trees.
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
module spamm_algebra_1d

  use spamm_real_precision
  use spamm_tree_1d
  use spamm_strings

#ifdef _OPENMP
  use omp_lib
#endif

#include "spamm_utility_macros.h"

  implicit none

  !> Interface for multiplication operations between different SpAMM types.
  !!
  !! The Sparse Approximate Matrix-Multiply (SpAMM):
  !! @f$ C \leftarrow A \cdot B @f$.
  interface multiply
     module procedure spamm_multiply_1d_by_1d
  end interface multiply

  !> Interface for trace operations.
  interface trace
  end interface trace

  !> Interface for additions operations between different SpAMM types.
  interface add
  end interface add

  !> Interface for filter operations (thresholding of small matrix
  !! elements).
  interface filter
  end interface filter

  !> Interface for norm operations.
  interface norm
  end interface norm

  !> Interface for dot product operations.
  interface dot
  end interface dot

contains

  !> Multiply two 1-D trees, i.e. a vector-vector inner product.
  !!
  !! Form the inner produce of two 1-D trees: \f$ C \leftarrow \alpha
  !! A \dot B \f$. The product is formed in place if possible. The
  !! matrix dimensions of the input vectors \f$ A \f$ and \f$ B \f$
  !! have to match.
  !!
  !! @param A Vector A.
  !! @param B Vector B.
  !! @param alpha The factor \f$ \alpha \f$.
  !! @param tolerance The SpAMM tolerance.
  function spamm_multiply_1d_by_1d (A, B, alpha, tolerance) result (C)

    type(tree_1d), pointer :: C
    type(tree_1d), pointer, intent(inout) :: A
    type(tree_1d), pointer, intent(in) :: B
    real(spamm_kind), optional, intent(in) :: alpha
    real(spamm_kind), optional, intent(in) :: tolerance

    real(spamm_kind) :: alpha_
    real(spamm_kind) :: tolerance_

    if(present(alpha)) then
       alpha_ = alpha
    else
       alpha_ = 0
    endif

    if(present(tolerance)) then
       tolerance_ = tolerance
    else
       tolerance_ = 0
    endif

    if(.not. associated(A)) then
       error stop
    endif

    if(.not. associated(B)) then
       error stop
    endif

    C => null()

  end function spamm_multiply_1d_by_1d

end module spamm_algebra_1d
