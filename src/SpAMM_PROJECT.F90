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

  USE  SpAMM_DERIVED
  USE  SpAMM_GLOBALS
  USE  SpAMM_MNGMENT
  USE  SpAMM_ALGEBRA

  IMPLICIT NONE

  !===============================================================================
  !  GLOBALS
  !===============================================================================

  !> Interface for spectral projections.
  INTERFACE RemapSpectralBounds201
     MODULE PROCEDURE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree
  END INTERFACE

  !> Interface for TC2.
  INTERFACE SpAMM_TC2
     MODULE PROCEDURE SpAMM_Quadratic_Trace_Correcting_Purification
  END INTERFACE

CONTAINS
  !=================================================================
  ! SPAMM ROUTINES FOR SPECTRAL PROJECTION
  !=================================================================
  ! Quadratic trace correcting purification
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SUBROUTINE SpAMM_Quadratic_Trace_Correcting_Purification(P,P2,Ne,TrP)
    TYPE(QuTree),POINTER   :: P,P2
    TYPE(QuTree),POINTER   :: T1
    REAL(SpAMM_KIND)       :: Ne,TrP,TrP2
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    TInitial=SpAMM_Get_Time()
    TrP=Trace(P)
    CALL Multiply(P,P,P2)                          ! P^2 <- P.P
    TrP2=Trace(P2)
    IF(ABS(TrP2-Ne)<ABS(SpAMM_Two*TrP-TrP2-Ne))THEN
       T1=>P; P=>P2; P2=>T1                        ! P <- P^2
    ELSE
       CALL Add(P,P2,SpAMM_Two,-SpAMM_One)         ! P <- 2*P-P^2
    ENDIF
    P%Norms=Norm(P)
    P%Norms%FrobeniusNorm=SQRT(P%Norms%FrobeniusNorm)
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Quadratic_Trace_Correcting_Purification",29)
  END SUBROUTINE SpAMM_Quadratic_Trace_Correcting_Purification

  !=================================================================
  ! REMAPING SPECTRAL BOUNDS TO [0,1]
  !=================================================================
  SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree(A)
    TYPE(QuTree),POINTER  :: A
    REAL(SpAMM_KIND)     :: RQIMin,RQIMax,SpectralExtent
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_MULTIPLY_THRESHOLD   =1D-7 !SpAMM_PRODUCT_TOLERANCE
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_CONVERGENCE_THRESHOLD=1D-3 !100d0*SpAMM_RQI_MULTIPLY_THRESHOLD
    REAL(SpAMM_KIND),PARAMETER  :: SpAMM_RQI_EVAL_ERROR_ESTIMATE  =2D-2 !100d0*SpAMM_RQI_CONVERGENCE_THRESHOLD
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
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
    A%Norms=Norm(A)
    A%Norms%FrobeniusNorm=Sqrt(A%Norms%FrobeniusNorm)
    TTotal=SpAMM_Get_Time()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree",31)
  END SUBROUTINE SpAMM_Remap_Spectral_Bounds_To_Zero_And_One_QuTree
  !=================================================================
  ! SPAMM ROUTINES FOR SPECTRAL ESTIMATION (EXTREMAL EIGENVALUES)
  !=================================================================
  ! RQI for finding the eigen bounds (Min/Max E.V.s)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Spectral_Bounds_Estimated_by_RQI_QuTree(A,RQIMin,RQIMax, &
             SpAMM_RQI_MULTIPLY_THRESHOLD,SpAMM_RQI_CONVERGENCE_THRESHOLD)
    INTEGER              :: N,I,CG
    REAL(SpAMM_KIND)     :: SpAMM_RQI_MULTIPLY_THRESHOLD, SpAMM_RQI_CONVERGENCE_THRESHOLD
    INTEGER, PARAMETER   :: NCG=1000
    TYPE(QuTree),POINTER :: A
    TYPE(BiTree),POINTER :: x=>NULL(),g=>NULL(),h=>NULL(),Ax=>NULL(), &
                           Ah=>NULL(),xOld=>NULL(),gOld=>NULL(),hOld=>NULL(), &
                           gTmp=>NULL(),hTmp=>NULL()
    REAL(SpAMM_KIND) ::  beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega, &
         xx,hh,xh,hx,xAx,xAh,hAx,hAh,GSize,xnorm
    REAL(SpAMM_KIND) :: RQIMin,RQIMax,Lambda,Lambdamax,LambdaMin
    !
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
       !
       CALL Multiply(x,xnorm)
       x%Norms=Norm(x)
       x%Norms%FrobeniusNorm=SQRT(x%Norms%FrobeniusNorm)
       WRITE(*,*)' XNORM = ',XNORM
!       STOP
       !
       CALL Multiply(h,SpAMM_Zero)
       CALL Multiply(g,SpAMM_Zero)
       CALL Multiply(xOld,SpAMM_Zero)
       CALL Multiply(hOld,SpAMM_Zero)
       CALL Multiply(gOld,SpAMM_Zero)
       !
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
          !
          IF(SQRT(Dot(g,g)/ABS(Omega))<SpAMM_RQI_CONVERGENCE_THRESHOLD.AND.CG>16)EXIT
          !
          IF(CG>1.AND.MOD(CG,15).NE.0)THEN
             !             beta=MAX(SpAMM_Zero,Dot(g,g-gOld)/Dot(gOld,gOld))
             beta=MAX(SpAMM_Zero,Dot(g,g)/Dot(gOld,gOld))
          ELSE
             beta=SpAMM_Zero
          ENDIF
          !
          ! h=g+beta*hOld
          CALL Add(h,SpAMM_One,g,beta,hOld)
          h%Norms=Norm(h)
          h%Norms%FrobeniusNorm=SQRT(h%Norms%FrobeniusNorm)
          !
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
          !
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
          !
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
!!$    REAL(SpAMM_DOUBLE), DIMENSION(N) :: Vals1
!!$    REAL(SpAMM_DOUBLE), DIMENSION(N,N) :: Temp
!!$    REAL(SpAMM_DOUBLE), DIMENSION(3*N) :: Work
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
!!$    REAL(SpAMM_DOUBLE), DIMENSION(N) :: Vals1
!!$    REAL(SpAMM_DOUBLE), DIMENSION(N,N) :: Temp
!!$    REAL(SpAMM_DOUBLE), DIMENSION(3*N) :: Work
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
