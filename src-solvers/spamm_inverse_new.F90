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
module spamm_inverse_new

#include "spamm_utility_macros.h"

  use spamm_tree_2d

  implicit none

contains
  !> Compute the inverse square root factorization of a matrix.
  !!
  !! @param S The matrix.
  !! @param Y The square root of the matrix, @f$ Y \leftarrow S^{1/2} @f$.
  !! @param Z The inverse square root of the matrix, @f$ Z \leftarrow S^{-1/2} @f$.
  !! @param tolerance The SpAMM tolerance.
  !! @param schulz_threshold The threshold for convergence. The error is
  !! defined as the @f$ \Vert X_{k} - I \Vert_{F} / (M N) @f$.
  subroutine spamm_inverse_sqrt_schulz_new (S, Y, Z, tolerance, schulz_threshold)

    use spamm_utilities

    type(tree_2d_symmetric), pointer, intent(in) :: S
    type(tree_2d_symmetric), pointer, intent(out) :: Y, Z
    real(spamm_kind), intent(in), optional :: tolerance
    real(spamm_kind), intent(in), optional :: schulz_threshold

    integer, parameter :: SCHULZ_ORDER = 3
    integer, parameter :: MAX_ITERATIONS = 200

    integer :: iteration

    ! Z_0 <- I
    !Z_0 => spamm_identity_2d(S)

    ! Y_0 <- S

    do iteration = 1, MAX_ITERATIONS
       ! X_{k} <- \lambda * Z^t_k * S * Z_k
       ! Error <- ||X_{k} - I||_{F}
       select case (SCHULZ_ORDER)
       case (2)
          ! Second order: T_{k} <- 1/2*(3*I-X_{k})
       case (3)
          ! Third order: T_{k} = 1/8*(15*I-10*X_{k}+3*X_{k}^{2})
       case (4)
          ! Fourth order: T_{k} = 1/16*(35*I-35*X_{k}+21*X_{k}^2-5*X({k}^3)
       case default
          LOG_FATAL("Schulz order "//to_string(SCHULZ_ORDER)//" is not implemented")
          error stop
       end select

    enddo

    ! Z_{k+1} <- Z_{k}*T_{k}
    ! Y_{k+1} <- T_{k}*Y_{k}

  end subroutine spamm_inverse_sqrt_schulz_new

end module spamm_inverse_new
