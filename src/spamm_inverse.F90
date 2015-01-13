!> Compute matrix inverses.
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
module spamm_inverse

  use spamm_real_precision

#include "spamm_utility_macros.h"

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
  subroutine spamm_inverse_sqrt_schulz (S, Y, Z, tolerance, schulz_threshold)

    use spamm_algebra
    use spamm_convert
    use spamm_management
    use spamm_types
    use spamm_utilities

    type(spamm_matrix_order_2), pointer, intent(in) :: S
    type(spamm_matrix_order_2), pointer, intent(inout) :: Y, Z
    real(spamm_kind), intent(in), optional :: tolerance
    real(spamm_kind), intent(in), optional :: schulz_threshold
    type(spamm_matrix_order_2), pointer :: X => null(), T => null(), &
      temp => null(), temp2 => null()

    integer, parameter :: SCHULZ_ORDER = 3
    integer, parameter :: MAX_ITERATIONS = 200
    integer :: iteration
    real(spamm_kind) :: lambda, error, local_tolerance, &
      local_schulz_threshold, number_operations


    if(present(tolerance)) then
      local_tolerance = tolerance
    else
      local_tolerance = 0
    endif

    if(present(schulz_threshold)) then
      local_schulz_threshold = schulz_threshold
    else
      local_schulz_threshold = 1e-8_spamm_kind
    endif

    if(associated(Y)) then
      call delete(Y)
    endif

    if(associated(Z)) then
      call delete(Z)
    endif

    LOG_INFO("Schulz -> approximate inverse sqrt")
    LOG_INFO("  tolerance = "//to_string(local_tolerance))
    LOG_INFO("  schulz threshold = "//to_string(local_schulz_threshold))

    ! Z_0 <- I
    Z => spamm_identity_matrix(S%M, S%N)

    ! Y_0 <- S
    call copy(S, Y)

    ! Something close to ideal scaling.
    lambda = 1/S%norm
    LOG_DEBUG("||S||_{F} = "//to_string(S%norm)//", lambda = "//to_string(lambda))

    ! Initialize error.
    error = 1

    ! Reset counters.
    number_operations = 0

    do iteration = 1, MAX_ITERATIONS
      ! X_{k} <- \lambda * Y_{k} * Z_{k}
      call multiply(Y, Z, X, local_tolerance)
      number_operations = number_operations+X%number_operations
      call multiply(X, lambda)
      number_operations = number_operations+X%number_operations

      ! Error <- ||X_{k} - I||_{F}
      T => spamm_identity_matrix(S%M, S%N)
      call add(T, X, -1.0_spamm_kind, 1.0_spamm_kind)
      error = T%norm/(T%M*T%N)
      call delete(T)

      LOG_INFO(to_string(iteration)//": error = "//to_string(error))

      if(error < local_schulz_threshold) then
        LOG_INFO("error below "//to_string(local_schulz_threshold))
        exit
      endif

      select case (SCHULZ_ORDER)
        case (2)
          ! Second order: T_{k} <- 1/2*(3*I-X_{k})
          T => spamm_identity_matrix(S%M, S%N)
          call add(T, X, 1.5_spamm_kind, -0.5_spamm_kind)

        case (3)
          ! Third order: T_{k} = 1/8*(15*I-10*X_{k}+3*X_{k}^{2})
          T => spamm_identity_matrix(S%M, S%N)
          call add(T, X, 15.0_spamm_kind, -10.0_spamm_kind)
          call multiply(X, X, temp, local_tolerance)
          call add(T, temp, 1.0_spamm_kind, 3.0_spamm_kind)
          call multiply(T, 1/8.0_spamm_kind)
          call delete(temp)

        case (4)
          ! Fourth order: T_{k} = 1/16*(35*I-35*X_{k}+21*X_{k}^2-5*X({k}^3)
          T => spamm_identity_matrix(S%M, S%N)
          call add(T, X, 35.0_spamm_kind, -35.0_spamm_kind)
          call multiply(X, X, temp, local_tolerance)
          call add(T, temp, 1.0_spamm_kind, 21.0_spamm_kind)
          call multiply(X, temp, temp2, local_tolerance)
          call add(T, temp2, 1.0_spamm_kind, -5.0_spamm_kind)
          call multiply(T, 1/16.0_spamm_kind)
          call delete(temp)
          call delete(temp2)

        case default
          LOG_FATAL("Schulz order "//to_string(SCHULZ_ORDER)//" is not implemented")
          error stop
      end select

      ! Z_{k+1} <- Z_{k}*T_{k}
      call multiply(Z, T, temp, local_tolerance)
      number_operations = number_operations+temp%number_operations
      call copy(temp, Z)
      call delete(temp)

      ! Y_{k+1} <- T_{k}*Y_{k}
      call multiply(T, Y, temp, local_tolerance)
      number_operations = number_operations+temp%number_operations
      call copy(temp, Y)
      call delete(temp)

      ! Free up T.
      call delete(T)
    enddo

    ! S^{-1/2} <- \sqrt{\lambda} Z_{k}
    call multiply(Z, sqrt(lambda))

    ! S^{1/2} <- \sqrt{\lambda} Y_{k}
    call multiply(Y, sqrt(lambda))

    LOG_INFO("Z (S^{-1/2} nnonzero = "//to_string(Z%number_nonzeros))
    LOG_INFO("  % fillin           = "//to_string(Z%number_nonzeros/(Z%M*Z%N)))
    LOG_INFO("Y (S^{+1/2} nnonzero = "//to_string(Y%number_nonzeros))
    LOG_INFO("  % fillin           = "//to_string(Y%number_nonzeros/(Y%M*Y%N)))
    LOG_INFO("number operations    = "//to_string(number_operations))

    call delete(X)
    call delete(T)

  end subroutine spamm_inverse_sqrt_schulz

end module spamm_inverse
