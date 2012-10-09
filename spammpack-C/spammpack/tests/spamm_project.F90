!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 3 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!-------------------------------------------------------------------------------
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SpAMMPack)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------

!> Projections of quadtree objects.
MODULE SpAMM_PROJECT

  USE spammpack

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

    TYPE(SpAMM_Matrix), INTENT(INOUT) :: P
    TYPE(SpAMM_Matrix), INTENT(INOUT) :: P2
    REAL(SpAMM_KIND), INTENT(IN)      :: Ne
    REAL(SpAMM_KIND), INTENT(OUT)     :: TrP

    TYPE(SpAMM_Matrix) :: T1
    REAL(SpAMM_KIND)   :: TrP2
    REAL(SpAMM_DOUBLE) :: TInitial, TTotal

    TInitial=SpAMM_Get_Time()

    TrP=Trace(P)
    CALL Multiply(P,P,P2)                          ! P^2 <- P.P
    TrP2=Trace(P2)

    IF(ABS(TrP2-Ne)<ABS(SpAMM_Two*TrP-TrP2-Ne))THEN
      T1=>P; P=>P2; P2=>T1                        ! P <- P^2
    ELSE
      CALL Add(P,P2,SpAMM_Two,-SpAMM_One)         ! P <- 2*P-P^2
    ENDIF

    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Quadratic_Trace_Correcting_Purification",29)

  END SUBROUTINE SpAMM_Quadratic_Trace_Correcting_Purification

  !> Remaping Spectral Bounds to [0, 1]
  !>
  !> @param A The matrix.
  SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree(A)

    TYPE(SpAMM_Matrix), intent(inout) :: A
    REAL(SpAMM_KIND)                  :: RQIMin,RQIMax,SpectralExtent
    REAL(SpAMM_KIND),PARAMETER        :: SpAMM_RQI_MULTIPLY_THRESHOLD   =1D-7 !SpAMM_PRODUCT_TOLERANCE
    REAL(SpAMM_KIND),PARAMETER        :: SpAMM_RQI_CONVERGENCE_THRESHOLD=1D-3 !100d0*SpAMM_RQI_MULTIPLY_THRESHOLD
    REAL(SpAMM_KIND),PARAMETER        :: SpAMM_RQI_EVAL_ERROR_ESTIMATE  =2D-2 !100d0*SpAMM_RQI_CONVERGENCE_THRESHOLD
    REAL(SpAMM_DOUBLE)                :: TInitial, TTotal

    TInitial=SpAMM_Get_Time()
    ! Find extremal eigenvalues by RQI
    CALL SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(A,RQIMin,RQIMax, &
      SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)
    ! Assume a 1% error in the estimated bounds
    RQIMin=RQIMin-SpAMM_RQI_EVAL_ERROR_ESTIMATE*ABS(RQIMin)
    RQIMax=RQIMax+SpAMM_RQI_EVAL_ERROR_ESTIMATE*ABS(RQIMax)
    SpectralExtent=RQIMax-RQIMin
    CALL Add(A,-RQIMax)
    CALL Multiply(A,-SpAMM_One/SpectralExtent)
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree",31)

  END SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree

  !> RQI for finding the eigen bounds (Min/Max E.V.s)
  !>
  !> @param A The matrix.
  SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(A,RQIMin,RQIMax, &
      SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)

    INTEGER                           :: I,CG
    REAL(SpAMM_KIND)                  :: SpAMM_RQI_MULTIPLY_THRESHOLD, SpAMM_RQI_CONVERGENCE_THRESHOLD
    INTEGER, PARAMETER                :: NCG=1000
    TYPE(SpAMM_Matrix), intent(inout) :: A
    TYPE(BiTree), POINTER             :: x=>NULL(),g=>NULL(),h=>NULL(),Ax=>NULL(),Ah=>NULL(),xOld=>NULL(),gOld=>NULL(),hOld=>NULL()
    REAL(SpAMM_KIND)                  :: beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega,xx,hh,xh,hx,xAx,xAh,hAx,hAh,xnorm
    REAL(SpAMM_KIND)                  :: RQIMin,RQIMax

    CALL New(x)
    CALL New(g)
    CALL New(h)
    CALL New(Ax)
    CALL New(Ah)
    CALL New(xOld)
    CALL New(gOld)
    CALL New(hOld)

    DO I=1,2
      IF(I==1)THEN
        CALL Copy(A,1,x)
      ELSE
        CALL Copy(A,SpAMM_MATRIX_DIMENSION,x)
      ENDIF
      xnorm=SpAMM_One/Dot(x,x)

      CALL Multiply(x,xnorm)
      x%Norms=Norm(x)
      x%Norms%FrobeniusNorm=SQRT(x%Norms%FrobeniusNorm)
      WRITE(*,*)' XNORM = ',XNORM
      ! STOP

      CALL Multiply(h,SpAMM_Zero)
      CALL Multiply(g,SpAMM_Zero)
      CALL Multiply(xOld,SpAMM_Zero)
      CALL Multiply(hOld,SpAMM_Zero)
      CALL Multiply(gOld,SpAMM_Zero)

      DO CG=1,NCG
        ! Intermediates
        xx=Dot(x,x)

        WRITE(*,*)xx
        WRITE(*,*)" NormA = ",Norm(A)

        CALL Multiply(A,x,Ax,SpAMM_RQI_MULTIPLY_THRESHOLD)
        xAx=Dot(x,Ax)
        omega=xAx/xx
        IF(I==1)THEN
          RQIMin=omega
        ELSE
          RQIMax=omega
        ENDIF

        WRITE(*,*)" omega = ",omega
        STOP

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
          !             beta=MAX(SpAMM_Zero,Dot(g,g-gOld)/Dot(gOld,gOld))
          beta=MAX(SpAMM_Zero,Dot(g,g)/Dot(gOld,gOld))
        ELSE
          beta=SpAMM_Zero
        ENDIF

        ! h=g+beta*hOld
        CALL Add(h,SpAMM_One,g,beta,hOld)
        h%Norms=Norm(h)
        h%Norms%FrobeniusNorm=SQRT(h%Norms%FrobeniusNorm)

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
        x%Norms=Norm(x)
        x%Norms%FrobeniusNorm=SQRT(x%Norms%FrobeniusNorm)

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

END MODULE SpAMM_PROJECT
