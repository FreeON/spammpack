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

  implicit none

contains

  !> Compute the inverse square root factorization of a matrix.
  !!
  !! @param S The matrix.
  !! @param tolerance The SpAMM tolerance.
  !!
  !! @return The inverse square root of the matrix, @f$ Z \leftarrow S^{-1/2} @f$.
  function spamm_inverse_sqrt_schulz (S, tolerance) result (Z)

    use spamm_algebra
    use spamm_convert
    use spamm_management
    use spamm_types
    use spamm_utilities

    type(spamm_matrix_2nd_order), pointer, intent(in) :: S
    real(spamm_kind), intent(in), optional :: tolerance
    type(spamm_matrix_2nd_order), pointer :: Z
    type(spamm_matrix_2nd_order), pointer :: X => null(), Y => null(), T => null(), temp => null()

    integer, parameter :: max_iterations = 20
    real(spamm_kind), parameter :: SCHULZ_THRESHOLD = 1e-10_spamm_kind
    integer :: iteration
    real(spamm_kind) :: lambda, error

    LOG_DEBUG("Schulz -> approximate inverse sqrt")

    ! Z_0 <- I
    Z => spamm_identity_matrix(S%M, S%N)

    ! Y_0 <- S
    call copy(S, Y)

    ! Something close to ideal scaling.
    lambda = 1/S%norm
    LOG_DEBUG("||S||_{F} = "//to_string(S%norm)//", lambda = "//to_string(lambda))

    ! Initialize error.
    error = 1

    do iteration = 1, max_iterations
      call print_spamm_2nd_order(Z, "Z_"//to_string(iteration))
      call print_spamm_2nd_order(Y, "Y_"//to_string(iteration))

      ! X_{k} <- \lambda * Y_{k} * Z_{k}
      call multiply(Y, Z, X, tolerance)
      call print_spamm_2nd_order(X, "X_"//to_string(iteration))
      call multiply(X, lambda)
      call print_spamm_2nd_order(X, "X_"//to_string(iteration))

      ! Error <- ||X_{k} - I||_{F}
      T => spamm_identity_matrix(S%M, S%N)
      call add(T, X, -1.0_spamm_kind, 1.0_spamm_kind)
      error = T%norm
      call delete(T)

      LOG_DEBUG(to_string(iteration)//": error = "//to_string(error))

      if(error < SCHULZ_THRESHOLD) then
        LOG_DEBUG("error ratio below "//to_string(SCHULZ_THRESHOLD))
        exit
      endif

      ! Second order: T_{k} <- 1/2*(3*I-X_{k})
      T => spamm_identity_matrix(S%M, S%N)
      call print_spamm_2nd_order(T, "T_"//to_string(iteration))
      call add(T, X, 3.0_spamm_kind, -1.0_spamm_kind)
      call print_spamm_2nd_order(T, "T_"//to_string(iteration))
      call multiply(T, 0.5_spamm_kind)
      call print_spamm_2nd_order(T, "T_"//to_string(iteration))

      ! Z_{k+1} <- Z_{k}*T_{k}
      call copy(Z, temp)
      call multiply(temp, T, Z, tolerance)

      ! Y_{k+1} <- T_{k}*Y_{k}
      call copy(Y, temp)
      call multiply(T, temp, Y, tolerance)

      ! Free up T.
      call delete(T)
    enddo

    ! S^{-1/2} <- \sqrt{\lambda} Z_{k}
    call multiply(Z, sqrt(lambda))

    call delete(X)
    call delete(Y)
    call delete(T)

  end function spamm_inverse_sqrt_schulz

end module spamm_inverse
