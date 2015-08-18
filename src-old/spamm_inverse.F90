!> Compute matrix inverses.
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
module spamm_inverse

  use spamm_algebra
  use spamm_convert
  use spamm_management
  use spamm_real_precision
  use spamm_types
  use spamm_utilities

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

    type(spamm_matrix_order_2), pointer, intent(in) :: S
    type(spamm_matrix_order_2), pointer, intent(inout) :: Y, Z
    real(spamm_kind), intent(in), optional :: tolerance
    real(spamm_kind), intent(in), optional :: schulz_threshold
    type(spamm_matrix_order_2), pointer :: X => null(), T => null(), &
      temp => null(), temp2 => null()

    integer, parameter :: SCHULZ_ORDER = 3
    integer, parameter :: MAX_ITERATIONS = 200
    integer :: iteration
    real(spamm_kind) :: lambda, error, error_last, rel_change, local_tolerance, &
      local_schulz_threshold, number_operations


    REAL(SpAMM_KIND) :: EvMin, EvMax
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_MULTIPLY_THRESHOLD   =1D-7 !SpAMM_PRODUCT_TOLERANCE
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_CONVERGENCE_THRESHOLD=1D-3 !100d0*SpAMM_RQI_MULTIPLY_THRESHOLD
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_EVAL_ERROR_ESTIMATE  =2D-2 !100d0*SpAMM_RQI_CONVERGENCE_THRESHOLD

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

    WRITE(*,*)' FINDING COND NUMBER = '
    ! Find extremal eigenvalues by RQI
    CALL SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(S%Root,EvMin,EvMax, &
      SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)

    WRITE(*,*)' MINMAX =',EvMin, EvMax
    WRITE(*,*)' COND S = ',EvMax/EvMin

    LOG_INFO("Schulz -> approximate inverse sqrt")
    LOG_INFO("  tolerance = "//to_string(local_tolerance))
    LOG_INFO("  schulz threshold = "//to_string(local_schulz_threshold))

    ! Z_0 <- I
    Z => spamm_identity_matrix(S%M, S%N)

    ! Y_0 <- S
    call copy(S, Y)

    ! Something close to ideal scaling.

    lambda = 2d0/(EvMin+EvMax)

    LOG_DEBUG("||S||_{F} = "//to_string(S%norm)//", lambda = "//to_string(lambda))

    ! Initialize error.
    error_last = 1


    ! Reset counters.
    number_operations = 0

    do iteration = 1, 17 !MAX_ITERATIONS
      ! X_{k} <- \lambda * Z^t_k * S * Z_k
      call multiply(S, Z, T, local_tolerance)
      number_operations = number_operations+T%number_operations
      call multiply( Z, T, X, local_tolerance)
      number_operations = number_operations+X%number_operations
      call multiply(X, lambda)
      number_operations = number_operations+X%number_operations
      ! Error <- ||X_{k} - I||_{F}
      T => spamm_identity_matrix(S%M, S%N)
      call add(T, X, -1.0_spamm_kind, 1.0_spamm_kind)
      error = T%norm/(T%M*T%N)
      rel_change=(error_last-error)/error
      error_last=error
      call delete(T)  !<<<<<<<<<<< probably overkill here
      WRITE(*,*)to_string(iteration)//"::: error = "//to_string(error)
      LOG_INFO(to_string(iteration)//": error = "//to_string(error))
      if( iteration > 5 .and. rel_change < 1d-2 ) then
        LOG_INFO("111 error below "//to_string(local_schulz_threshold))
        !        exit
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

!      ! Y_{k+1} <- T_{k}*Y_{k}
!      call multiply(T, Y, temp, local_tolerance)
!      number_operations = number_operations+temp%number_operations
!      call copy(temp, Y)
!      call delete(temp)

      ! Free up T.
      call delete(T)
    enddo


    goto 9999

!!$    ! Reset counters.
!!$    number_operations = 0
!!$    do iteration = 1, MAX_ITERATIONS
!!$      ! X_{k} <- \lambda * Y_{k} * Z_{k}
!!$      call multiply(Y, Z, X, local_tolerance)
!!$      number_operations = number_operations+X%number_operations
!!$      call multiply(X, lambda)
!!$      number_operations = number_operations+X%number_operations
!!$      ! Error <- ||X_{k} - I||_{F}
!!$      T => spamm_identity_matrix(S%M, S%N)
!!$      call add(T, X, -1.0_spamm_kind, 1.0_spamm_kind)
!!$      error = T%norm/(T%M*T%N)
!!$      rel_change=(error_last-error)/error
!!$      error_last=error
!!$      call delete(T)
!!$      WRITE(*,*)to_string(iteration)//": error = "//to_string(error)
!!$      LOG_INFO(to_string(iteration)//": error = "//to_string(error))
!!$      if( iteration > 5 .and. rel_change < 1d-2 ) then
!!$!      if(error < local_schulz_threshold) then
!!$        LOG_INFO("error below "//to_string(local_schulz_threshold))
!!$        exit
!!$      endif
!!$      select case (SCHULZ_ORDER)
!!$        case (2)
!!$          ! Second order: T_{k} <- 1/2*(3*I-X_{k})
!!$          T => spamm_identity_matrix(S%M, S%N)
!!$          call add(T, X, 1.5_spamm_kind, -0.5_spamm_kind)
!!$        case (3)
!!$          ! Third order: T_{k} = 1/8*(15*I-10*X_{k}+3*X_{k}^{2})
!!$          T => spamm_identity_matrix(S%M, S%N)
!!$          call add(T, X, 15.0_spamm_kind, -10.0_spamm_kind)
!!$          call multiply(X, X, temp, local_tolerance)
!!$          call add(T, temp, 1.0_spamm_kind, 3.0_spamm_kind)
!!$          call multiply(T, 1/8.0_spamm_kind)
!!$          call delete(temp)
!!$        case (4)
!!$          ! Fourth order: T_{k} = 1/16*(35*I-35*X_{k}+21*X_{k}^2-5*X({k}^3)
!!$          T => spamm_identity_matrix(S%M, S%N)
!!$          call add(T, X, 35.0_spamm_kind, -35.0_spamm_kind)
!!$          call multiply(X, X, temp, local_tolerance)
!!$          call add(T, temp, 1.0_spamm_kind, 21.0_spamm_kind)
!!$          call multiply(X, temp, temp2, local_tolerance)
!!$          call add(T, temp2, 1.0_spamm_kind, -5.0_spamm_kind)
!!$          call multiply(T, 1/16.0_spamm_kind)
!!$          call delete(temp)
!!$          call delete(temp2)
!!$        case default
!!$          LOG_FATAL("Schulz order "//to_string(SCHULZ_ORDER)//" is not implemented")
!!$          error stop
!!$      end select
!!$
!!$      ! Z_{k+1} <- Z_{k}*T_{k}
!!$      call multiply(Z, T, temp, local_tolerance)
!!$      number_operations = number_operations+temp%number_operations
!!$      call copy(temp, Z)
!!$      call delete(temp)
!!$
!!$      ! Y_{k+1} <- T_{k}*Y_{k}
!!$      call multiply(T, Y, temp, local_tolerance)
!!$      number_operations = number_operations+temp%number_operations
!!$      call copy(temp, Y)
!!$      call delete(temp)
!!$
!!$      ! Free up T.
!!$      call delete(T)
!!$    enddo
!!$

9999 continue

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


  !> Spamm routines for spectral estimation (extremal eigenvalues)
  !!
  !! RQI for finding the eigen bounds (Min/Max E.V.s).
  SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(A,RQIMin,RQIMax, &
      SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)

    integer              :: I,CG
    real(SpAMM_KIND)     :: SpAMM_RQI_MULTIPLY_THRESHOLD, SpAMM_RQI_CONVERGENCE_THRESHOLD
    integer, parameter   :: NCG=1000
    type(qutree), pointer, intent(in) :: A
    type(BiTree),pointer :: x=>NULL(),g=>NULL(),h=>NULL(),Ax=>NULL(), &
         Ah=>NULL(),xOld=>NULL(),gOld=>NULL(),hOld=>NULL()
    real(SpAMM_KIND)     :: beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega, &
         xx,hh,xh,hx,xAx,xAh,hAx,hAh,xnorm
    real(SpAMM_KIND)     :: RQIMin, RQIMax

    CALL New(x, A%i_lower, A%i_upper)
    CALL New(g, A%i_lower, A%i_upper)
    CALL New(h, A%i_lower, A%i_upper)
    CALL New(Ax, A%i_lower, A%i_upper)
    CALL New(Ah, A%i_lower, A%i_upper)
    CALL New(xOld, A%i_lower, A%i_upper)
    CALL New(gOld, A%i_lower, A%i_upper)
    CALL New(hOld, A%i_lower, A%i_upper)

    DO I=1,2
!      IF(I==1)THEN
        CALL Copy(A,1,x)
!      ELSE
!         WRITE(*,*)' I=2 copy ....'

!        CALL Copy(A,1,x) !SpAMM_MATRIX_DIMENSION,x)
!      ENDIF

      write(*,*)' 3'

      xnorm=SpAMM_One/Dot(x,x)
      CALL Multiply(x,xnorm)
      ! This call should be redundant  <<<<<<<<<<<>>>>>>>>>>>>>>
      x%Norm=SQRT(Norm(x))
      CALL Multiply(h, 0d0)
      CALL Multiply(g, 0d0)
      CALL Multiply(xOld, 0d0)
      CALL Multiply(hOld, 0d0)
      CALL Multiply(gOld, 0d0)
      DO CG=1,NCG
        ! Intermediates
        xx=Dot(x,x)
        CALL Multiply(A,x,Ax,SpAMM_RQI_MULTIPLY_THRESHOLD)
        xAx=Dot(x,Ax)
        omega=xAx/xx
        IF(I==1)THEN
          RQIMin=omega
        ELSE
          RQIMax=omega
        ENDIF
        ! Gradient of extremal quotients (eigenvalues): one is + the other is - ...
        IF(I==1)THEN
          ! g=2*(Ax-omega*x)/xx
          CALL Add(g,SpAMM_Two/xx,Ax,-SpAMM_Two*omega/xx,x)
        ELSE
          ! g=-2*(Ax-omega*x)/xx
          CALL Add(g,-SpAMM_Two/xx,Ax,SpAMM_Two*omega/xx,x)
        ENDIF

        IF(SQRT(Dot(g,g)/ABS(Omega))<SpAMM_RQI_CONVERGENCE_THRESHOLD.AND.CG>16)EXIT

        IF(CG>1.AND.MOD(CG,15).NE.0)THEN
          !             beta=MAX(0,Dot(g,g-gOld)/Dot(gOld,gOld))
          beta=MAX(0d0, Dot(g,g)/Dot(gOld,gOld))
        ELSE
          beta = 0
        ENDIF

        ! h=g+beta*hOld
        CALL Add(h,SpAMM_One,g,beta,hOld)
        h%Norm=SQRT(Norm(h))
        CALL Copy(g,gOld)
        CALL Copy(h,hOld)
        ! Ah=A.h
        CALL Multiply(A,h,Ah,SpAMM_RQI_MULTIPLY_THRESHOLD)
        hx =Dot(h,x)
        hh =Dot(h,h)
        xAh=Dot(x,Ah)
        hAx=Dot(h,Ax)
        hAh=Dot(h,Ah)
        ! By symmetry
        xh=hx
        hAx=xAh
        LambdaPlus=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx+SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
          -SpAMM_Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
          /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
        LambdaMins=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx-SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
          -SpAMM_Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
          /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
        !
        RQIPlus=(xAx+LambdaPlus*(xAh+hAx)+hAh*LambdaPlus**2) &
          /( xx+LambdaPlus*(xh+hx)  +hh *LambdaPlus**2)
        RQIMins=(xAx+LambdaMins*(xAh+hAx)+hAh*LambdaMins**2) &
          /( xx+LambdaMins*(xh+hx)  +hh *LambdaMins**2)

        IF(I==1)THEN
          IF(RQIMins<RQIPlus)THEN
            ! x=x+LambdaMins*h
            CALL Add(x,SpAMM_One,h,LambdaMins)
          ELSE
            ! x=x+LambdaPlus*h
            CALL Add(x,SpAMM_One,h,LambdaPlus)
          ENDIF
        ELSE
          IF(RQIMins>RQIPlus)THEN
            ! x=x+LambdaMins*h
            CALL Add(x,SpAMM_One,h,LambdaMins)
          ELSE
            ! x=x+LambdaPlus*h
            CALL Add(x,SpAMM_One,h,LambdaPlus)
          ENDIF
        ENDIF
        x%Norm=SQRT(Norm(x))
        IF(I==1)THEN
          WRITE(*,33)omega,SQRT(Dot(g,g))/ABS(Omega),CG
        ELSE
          WRITE(*,44)omega,SQRT(Dot(g,g))/ABS(Omega),CG
        ENDIF
      END DO
      IF(I==1)THEN
        WRITE(*,33)omega,SQRT(Dot(g,g))/ABS(Omega),CG
      ELSE
        WRITE(*,44)omega,SQRT(Dot(g,g))/ABS(Omega),CG
      ENDIF
    ENDDO

33  FORMAT(' MIN E.V. = ',E10.3,', GRAD RQI = ',E10.3,' in ',I4,' NLCG steps')
44  FORMAT(' MAX E.V. = ',E10.3,', GRAD RQI = ',E10.3,' in ',I4,' NLCG steps')

    CALL Delete(x)
    CALL Delete(g)
    CALL Delete(h)
    CALL Delete(Ax)
    CALL Delete(Ah)
    CALL Delete(xOld)
    CALL Delete(gOld)
    CALL Delete(hOld)

  END SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree


end module spamm_inverse
