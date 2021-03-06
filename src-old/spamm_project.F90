!> Projections of quadtree objects.
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
module SpAMM_PROJECT

  use spamm_algebra
  use spamm_globals
  use spamm_management
  use spamm_real_precision
  use spamm_types
  use spamm_utilities

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RemapSpectralBounds201
  PUBLIC :: SpAMM_TC2

  !> Interface for spectral projections.
  INTERFACE RemapSpectralBounds201
    MODULE PROCEDURE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree
  END INTERFACE

  !> Interface for TC2.
  INTERFACE SpAMM_TC2
    MODULE PROCEDURE SpAMM_Quadratic_Trace_Correcting_Purification
  END INTERFACE

CONTAINS

  !> Quadratic trace correcting purification.
  !!
  !! @param P The density matrix.
  !! @param P2 The squared density matrix.
  !! @param Ne The number of electrons.
  !! @param TrP The trace of P.
  SUBROUTINE SpAMM_Quadratic_Trace_Correcting_Purification (P, P2, Ne, TrP)
    TYPE(QuTree), POINTER, INTENT(INOUT) :: P
    TYPE(QuTree), POINTER, INTENT(INOUT) :: P2
    REAL(SpAMM_KIND), INTENT(IN)         :: Ne
    REAL(SpAMM_KIND), INTENT(OUT)        :: TrP
    TYPE(QuTree), POINTER :: T1
    REAL(SpAMM_KIND)      :: TrP2
    REAL(kind(0d0))    :: TInitial, TTotal

    TInitial=SpAMM_Get_Time()
    TrP=Trace(P)
    CALL Multiply(P,P,P2)                          ! P^2 <- P.P
    TrP2=Trace(P2)
    IF(ABS(TrP2-Ne)<ABS(SpAMM_Two*TrP-TrP2-Ne))THEN
      T1=>P; P=>P2; P2=>T1                        ! P <- P^2
    ELSE
      CALL Add(P,P2,SpAMM_Two,-SpAMM_One)         ! P <- 2*P-P^2
    ENDIF
    P%Norm=SQRT(Norm(P))
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Quadratic_Trace_Correcting_Purification",29)
  END SUBROUTINE SpAMM_Quadratic_Trace_Correcting_Purification

  !> Remaping spectral bounds to \f$ [0, 1] \f$.
  !!
  !! @bug Commented out "call add(A, -RQMax)" statement.
  !!
  !! @param A Pointer to matrix.
  SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree (A)

    TYPE(qutree), POINTER, intent(inout) :: A
    REAL(SpAMM_KIND) :: RQIMin,RQIMax,SpectralExtent
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_MULTIPLY_THRESHOLD   =1D-7 !SpAMM_PRODUCT_TOLERANCE
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_CONVERGENCE_THRESHOLD=1D-3 !100d0*SpAMM_RQI_MULTIPLY_THRESHOLD
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_EVAL_ERROR_ESTIMATE  =2D-2 !100d0*SpAMM_RQI_CONVERGENCE_THRESHOLD
    REAL(Kind(0d0)) :: TInitial, TTotal

    TInitial=SpAMM_Get_Time()
    ! Find extremal eigenvalues by RQI
    CALL SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(A,RQIMin,RQIMax, &
      SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)
    ! Assume a 1% error in the estimated bounds
    RQIMin=RQIMin-SpAMM_RQI_EVAL_ERROR_ESTIMATE*ABS(RQIMin)
    RQIMax=RQIMax+SpAMM_RQI_EVAL_ERROR_ESTIMATE*ABS(RQIMax)
    SpectralExtent=RQIMax-RQIMin
    !CALL Add(A, -RQIMax)
    CALL Multiply(A, -SpAMM_One/SpectralExtent)
    A%Norm = SQRT(Norm(A))
    TTotal = SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree",31)

  END SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree

  !> Spamm routines for spectral estimation (extremal eigenvalues)
  !!
  !! RQI for finding the eigen bounds (Min/Max E.V.s).
  SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(A,RQIMin,RQIMax, &
      SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)

    INTEGER              :: I,CG
    REAL(SpAMM_KIND)     :: SpAMM_RQI_MULTIPLY_THRESHOLD, SpAMM_RQI_CONVERGENCE_THRESHOLD
    INTEGER, PARAMETER   :: NCG=1000
    TYPE(qutree), POINTER, intent(in) :: A
    TYPE(BiTree),POINTER :: x=>NULL(),g=>NULL(),h=>NULL(),Ax=>NULL(), &
      Ah=>NULL(),xOld=>NULL(),gOld=>NULL(),hOld=>NULL()
    REAL(SpAMM_KIND)     :: beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega, &
      xx,hh,xh,hx,xAx,xAh,hAx,hAh,xnorm
    REAL(SpAMM_KIND)     :: RQIMin,RQIMax

    CALL New(x, A%i_lower, A%i_upper)
    CALL New(g, A%i_lower, A%i_upper)
    CALL New(h, A%i_lower, A%i_upper)
    CALL New(Ax, A%i_lower, A%i_upper)
    CALL New(Ah, A%i_lower, A%i_upper)
    CALL New(xOld, A%i_lower, A%i_upper)
    CALL New(gOld, A%i_lower, A%i_upper)
    CALL New(hOld, A%i_lower, A%i_upper)

    DO I=1,2
      IF(I==1)THEN
        CALL Copy(A,1,x)
      ELSE
        CALL Copy(A,SpAMM_MATRIX_DIMENSION,x)
      ENDIF
      xnorm=SpAMM_One/Dot(x,x)
      CALL Multiply(x,xnorm)
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
          !             beta=MAX(0d0, Dot(g,g-gOld)/Dot(gOld,gOld))
          beta=MAX(0d0, Dot(g,g)/Dot(gOld,gOld))
        ELSE
          beta=0
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

end module SpAMM_PROJECT

!!$
!!$  SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_SpAMM_Zero_And_SpAMM_One(N,A)
!!$    INTEGER                              :: N
!!$    REAL(SpAMM_KIND),DIMENSION(1:N,1:N)  :: A
!!$    REAL(SpAMM_KIND) :: RQIMin,RQIMax,SpectralExtent
!!$    CALL SpAMM_Spectral_Bounds_Estimated_by_RQI(N,A,RQIMin,RQIMax)
!!$    ! Assume a 1% error in the estimated bounds
!!$    RQIMin=RQIMin-1D-2*ABS(RQIMin)
!!$    RQIMax=RQIMax+1D-2*ABS(RQIMax)
!!$    SpectralExtent=RQIMax-RQIMin
!!$!    RETURN
!!$    DO I=1,N
!!$       A(I,I)=A(I,I)-RQIMax
!!$    ENDDO
!!$    A=(-SpAMM_One/SpectralExtent)*A
!!$
!!$  END SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_SpAMM_Zero_And_SpAMM_One
!!$
!!$
!!$
!!$  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$  ! Spectral bounds through RQI
!!$  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$  SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI(N,A,RQIMin,RQIMax)
!!$    INTEGER :: N,I,CG,NCG
!!$    REAL(SpAMM_KIND),DIMENSION(1:N,1:N)  ::  A
!!$    REAL(SpAMM_KIND),DIMENSION(1:N)      ::  x,g,h,Ax,Ah,xOld,gOld,hOld
!!$    REAL(SpAMM_KIND)                     ::  beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega, &
!!$         xx,hh,xh,hx,xAx,xAh,hAx,hAh,GSize
!!$    REAL(SpAMM_KIND) :: RQIMin,RQIMax,Lambda,Lambdamax,LambdaMin
!!$
!!$    INTEGER :: LWork,Info
!!$    REAL(Kind(0d0)), DIMENSION(N) :: Vals1
!!$    REAL(Kind(0d0)), DIMENSION(N,N) :: Temp
!!$    REAL(Kind(0d0)), DIMENSION(3*N) :: Work
!!$    LWork=MAX(1,3*N)
!!$    Vals2=0.0D0
!!$    Temp=A
!!$    Call DSYEV('V','U',N,Temp,N,Vals1,Work,LWork,Info)
!!$    WRITE(*,*)' MIN/MAX EigenValues = ',Vals1(1),Vals1(N)
!!$
!!$    RQIMin=Vals1(1)
!!$    RQIMax=Vals1(N)
!!$    RETURN
!!$
!!$    NCG=50
!!$
!!$    !    X(:,1)=Temp(:,1)
!!$    !    X(:,2)=Temp(:,N)
!!$
!!$    DO I=1,2
!!$       IF(I==1)THEN
!!$          x=A(:,1)
!!$       ELSE
!!$          x=A(:,N)
!!$       ENDIF
!!$
!!$
!!$!       WRITE(*,*)' xnorm = ',DOT_PRODUCT(x,x)
!!$
!!$       x=x/DOT_PRODUCT(x,x)
!!$
!!$       h=SpAMM_Zero
!!$       g=SpAMM_Zero
!!$       xOld=SpAMM_Zero
!!$       hOld=SpAMM_Zero
!!$       gOld=SpAMM_Zero
!!$       DO CG=1,NCG
!!$          ! Intermediates
!!$          xx=DOT_PRODUCT(x,x)
!!$          Ax=MATMUL(A,x)
!!$          xAx=DOT_PRODUCT(x,Ax)
!!$
!!$
!!$          omega=xAx/xx
!!$          IF(I==1)THEN
!!$             RQIMin=omega
!!$          ELSE
!!$             RQIMax=omega
!!$          ENDIF
!!$          ! Gradient of extremal quotients (eigenvalues): one is + the other is - ...
!!$          IF(I==1)THEN
!!$             g= SpAMM_Two*(Ax-omega*x)/xx
!!$          ELSE
!!$             g=-SpAMM_Two*(Ax-omega*x)/xx
!!$          ENDIF
!!$          !
!!$          IF(DOT_PRODUCT(g,g)<1D-8.AND.CG>20)EXIT
!!$          !
!!$          IF(CG>1.AND.MOD(CG,15).NE.0)THEN
!!$!             beta=MAX(SpAMM_Zero,DOT_PRODUCT(g,g-gOld)/DOT_PRODUCT(gOld,gOld))
!!$             beta=MAX(SpAMM_Zero,DOT_PRODUCT(g,g)/DOT_PRODUCT(gOld,gOld))
!!$          ELSE
!!$             beta=SpAMM_Zero
!!$          ENDIF
!!$          !
!!$          h=g+beta*hOld
!!$
!!$!          WRITE(*,*)' H = ',dot_product(h,h)
!!$
!!$          hOld=h
!!$          gOld=g
!!$          !
!!$          Ah =MATMUL(A,h)
!!$          hx =DOT_PRODUCT(h,x)
!!$          hh =DOT_PRODUCT(h,h)
!!$          xAh=DOT_PRODUCT(x,Ah)
!!$          hAx=DOT_PRODUCT(h,Ax)
!!$          hAh=DOT_PRODUCT(h,Ah)
!!$
!!$!          WRITE(*,*)' xx = ',xx
!!$!          WRITE(*,*)' hh = ',hh
!!$!          WRITE(*,*)' hx = ',hx
!!$!          WRITE(*,*)' hAx = ',hAx
!!$!          WRITE(*,*)' hAh = ',hAh
!!$
!!$
!!$
!!$          !
!!$
!!$
!!$          xh=hx
!!$          hAx=xAh
!!$          !
!!$          LambdaPlus=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx+SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
!!$               -Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
!!$               /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
!!$          LambdaMins=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx-SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
!!$               -Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
!!$               /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
!!$          !
!!$          RQIPlus=(xAx+LambdaPlus*(xAh+hAx)+hAh*LambdaPlus**2) &
!!$               /( xx+LambdaPlus*(xh+hx)  +hh *LambdaPlus**2)
!!$          RQIMins=(xAx+LambdaMins*(xAh+hAx)+hAh*LambdaMins**2) &
!!$               /( xx+LambdaMins*(xh+hx)  +hh *LambdaMins**2)
!!$
!!$          IF(I==1)THEN
!!$             IF(RQIMins<RQIPlus)THEN
!!$                x=x+LambdaMins*h
!!$             ELSE
!!$                x=x+LambdaPlus*h
!!$             ENDIF
!!$          ELSE
!!$             IF(RQIMins>RQIPlus)THEN
!!$                x=x+LambdaMins*h
!!$             ELSE
!!$                x=x+LambdaPlus*h
!!$             ENDIF
!!$          ENDIF
!!$
!!$!          WRITE(*,*)' Lambda = ',LambdaPlus,LambdaMins
!!$
!!$!          WRITE(*,*)'AFTERX = ',DOT_PRODUCT(X,X)
!!$!          STOP
!!$
!!$
!!$       END DO
!!$       IF(I==1)THEN
!!$          !          WRITE(*,33)omega,vals1(1),CG
!!$          WRITE(*,33)omega,SpAMM_Zero,CG
!!$       ELSE
!!$          !          WRITE(*,44)omega,vals1(N),CG
!!$          WRITE(*,44)omega,SpAMM_Zero,CG
!!$       ENDIF
!!$    ENDDO
!!$33  FORMAT(' MIN EIGEN VALUE = ',E18.6,' TRUE VAL = ',E18.6,' in ',I3,' RQI steps')
!!$44  FORMAT(' MAX EIGEN VALUE = ',E18.6,' TRUE VAL = ',E18.6,' in ',I3,' RQI steps')
!!$  END SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI
!!$
!!$  SUBROUTINE SpAMM_RQI_min_max(N,A)
!!$    REAL(SpAMM_KIND),DIMENSION(1:N,1:N)  ::  A
!!$    REAL(SpAMM_KIND),DIMENSION(1:N)      ::  x,g,h,Ax,Ah,xOld,gOld,hOld
!!$    REAL(SpAMM_KIND)                     ::  beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega, &
!!$         xx,hh,xh,hx,xAx,xAh,hAx,hAh,GSize
!!$    INTEGER :: I,CG,NCG
!!$    REAL(SpAMM_KIND) :: RQIMin,RQIMax,Lambda,Lambdamax,LambdaMin
!!$
!!$    INTEGER :: LWork,Info
!!$    REAL(Kind(0d0)), DIMENSION(N) :: Vals1
!!$    REAL(Kind(0d0)), DIMENSION(N,N) :: Temp
!!$    REAL(Kind(0d0)), DIMENSION(3*N) :: Work
!!$
!!$    LOGICAL :: MinMaxTest
!!$
!!$    LWork=MAX(1,3*N)
!!$    Vals2=0.0D0
!!$    Temp=A
!!$    Call DSYEV('V','U',N,Temp,N,Vals1,Work,LWork,Info)
!!$    WRITE(*,*)' MIN/MAX EigenValues = ',Vals1(1),Vals1(N)
!!$
!!$    NCG=50
!!$
!!$    !    X(:,1)=Temp(:,1)
!!$    !    X(:,2)=Temp(:,N)
!!$
!!$    DO I=1,2
!!$       IF(I==1)THEN
!!$          x=A(:,1)
!!$       ELSE
!!$          x=A(:,N)
!!$       ENDIF
!!$       x=x/DOT_PRODUCT(x,x)
!!$       h=SpAMM_Zero
!!$       g=SpAMM_Zero
!!$       xOld=SpAMM_Zero
!!$       hOld=SpAMM_Zero
!!$       gOld=SpAMM_Zero
!!$       DO CG=1,NCG
!!$          ! Intermediates
!!$          xx=DOT_PRODUCT(x,x)
!!$          Ax=MATMUL(A,x)
!!$          xAx=DOT_PRODUCT(x,Ax)
!!$          omega=xAx/xx
!!$          ! Gradient of extremal quotients (eigenvalues): one is + the other is - ...
!!$          IF(I==1)THEN
!!$             g= SpAMM_Two*(Ax-omega*x)/xx
!!$          ELSE
!!$             g=-SpAMM_Two*(Ax-omega*x)/xx
!!$          ENDIF
!!$          !
!!$          IF(DOT_PRODUCT(g,g)<1D-8.AND.CG>20)EXIT
!!$          !
!!$          IF(CG>1.AND.MOD(CG,15).NE.0)THEN
!!$             beta=MAX(SpAMM_Zero,DOT_PRODUCT(g,g-gOld)/DOT_PRODUCT(gOld,gOld))
!!$          ELSE
!!$             beta=SpAMM_Zero
!!$          ENDIF
!!$          !
!!$          h=g+beta*hOld
!!$          hOld=h
!!$          gOld=g
!!$          !
!!$          Ah =MATMUL(A,h)
!!$          hx =DOT_PRODUCT(h,x)
!!$          hh =DOT_PRODUCT(h,h)
!!$          xAh=DOT_PRODUCT(x,Ah)
!!$          hAx=DOT_PRODUCT(h,Ax)
!!$          hAh=DOT_PRODUCT(h,Ah)
!!$          !
!!$          xh=hx
!!$          hAx=xAh
!!$          !
!!$          LambdaPlus=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx+SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
!!$               -Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
!!$               /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
!!$          LambdaMins=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx-SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
!!$               -Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
!!$               /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
!!$          !
!!$          RQIPlus=(xAx+LambdaPlus*(xAh+hAx)+hAh*LambdaPlus**2) &
!!$               /( xx+LambdaPlus*(xh+hx)  +hh *LambdaPlus**2)
!!$          RQIMins=(xAx+LambdaMins*(xAh+hAx)+hAh*LambdaMins**2) &
!!$               /( xx+LambdaMins*(xh+hx)  +hh *LambdaMins**2)
!!$
!!$          IF(I==1)THEN
!!$             IF(RQIMins<RQIPlus)THEN
!!$                x=x+LambdaMins*h
!!$             ELSE
!!$                x=x+LambdaPlus*h
!!$             ENDIF
!!$          ELSE
!!$             IF(RQIMins>RQIPlus)THEN
!!$                x=x+LambdaMins*h
!!$             ELSE
!!$                x=x+LambdaPlus*h
!!$             ENDIF
!!$          ENDIF
!!$
!!$       END DO
!!$       IF(I==1)THEN
!!$          WRITE(*,33)omega,vals1(1),CG
!!$       ELSE
!!$          WRITE(*,44)omega,vals1(N),CG
!!$       ENDIF
!!$    ENDDO
!!$33  FORMAT(' MIN EIGEN VALUE = ',E18.6,' TRUE VAL = ',E18.6,' in ',I3,' RQI steps')
!!$44  FORMAT(' MAX EIGEN VALUE = ',E18.6,' TRUE VAL = ',E18.6,' in ',I3,' RQI steps')
!!$  END SUBROUTINE SpAMM_RQI_min_max
