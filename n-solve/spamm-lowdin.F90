!! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan"

!!$module GSOLVE_ALGEBRA
!!$
!!$    use spamm_algebra
!!$    use spamm_convert
!!$    use spamm_management
!!$    use spamm_types
!!$    use spamm_utilities
!!$
!!$    implicit none
!!$    INTERFACE Id
!!$       MODULE PROCEDURE GSOLVE_ALGEBRA_Initialize_Identity
!!$       MODULE PROCEDURE GSOLVE_ALGEBRA_Allocated_Identity
!!$    END INTERFACE Id
!!$    
!!$contains
!!$ 
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$  !! IT should be implicit that this multiply is accumulative.  That is, is assuming
!!$  !! A <- A*B or B <- A*B.  This is all fun and games until you go OMP.
!!$!  FUNCTION GSOLVE_MULTIPLY_QuTree_x_QuTree(qA, qB, Tolerance, Transpose, qC) RESULT C
!!$!    TYPE(QuTree),               POINTER, INTENT(IN) :: qA, qB
!!$!    TYPE(QuTree),     OPTIONAL, POINTER, INTENT(IN) :: qC
!!$!    REAL(SpAMM_KIND), OPTIONAL                      :: Threshold
!!$!    LOGICAL,          OPTIONAL                      :: Transpose
!!$!
!!$!    CALL SpAMM_Multiply_QuTree_x_QuTree(qA, qB, qC, Threshold, Transpose)
!!$!  END FUNCTION GSOLVE_MULTIPLY
!!$
!!$
!!$
!!$END module GSOLVE_ALGEBRA

MODULE GSOLVE_RQI
  USE SpAMM_algebra
  USE SpAMM_convert
  USE SpAMM_management
  USE SpAMM_types
  USE SpAMM_utilities
  USE test_utilities
CONTAINS
  !> SpAMM routines for spectral estimation (extremal eigenvalues)
  !!
  !! RQI for finding the eigen bounds (Min/Max E.V.s).
  SUBROUTINE GSOLVE_RQI_Extrema_RQI(A,RQIMin,RQIMax,MinVec,MaxVec,Tau,CnvrgCrit)

    TYPE(spamm_matrix_order_2) , POINTER, INTENT(IN)    :: A
!    TYPE(QuTree), POINTER, INTENT(IN)    :: A
    REAL(SpAMM_KIND),      INTENT(OUT)   :: RQIMin,RQIMax
    REAL(SpAMM_KIND),      INTENT(IN)    :: Tau, CnvrgCrit
    TYPE(BiTree), POINTER, INTENT(INOUT) :: MinVec,MaxVec

    INTEGER              :: I,CG, MM
    INTEGER, PARAMETER   :: NCG=5000
    TYPE(BiTree), POINTER :: x=>NULL(), g=>NULL(), h=>NULL(), Ax=>NULL(), &
                             Ah=>NULL(),xOld=>NULL(),gOld=>NULL(),hOld=>NULL()
  
    REAL(SpAMM_KIND)     :: beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega,omega_old
    REAL(SpAMM_KIND)     :: xx,hh,xh,hx,xAx,xAh,hAx,hAh,xnorm,dot_old

    !--------------------------------------------------------------------------------------
    real(spamm_kind), dimension(:, :), allocatable :: A_dense
    real(spamm_kind), dimension(:), allocatable :: MinVec_dense,MaxVec_dense, &
                       x_dense, g_dense, h_dense, Ax_dense,Ah_dense,xOld_dense,gOld_dense,hOld_dense

    REAL(SpAMM_KIND)     :: beta_dense,LambdaPlus_dense,LambdaMins_dense,RQIPlus_dense,RQIMins_dense,omega_dense
    REAL(SpAMM_KIND)     :: xx_dense,hh_dense,xh_dense,hx_dense,xAx_dense,xAh_dense,hAx_dense,hAh_dense,xnorm_dense
    REAL(SpAMM_KIND)     :: RQIMin_dense,RQIMax_dense,dot_old_dense

    MM=A%Root%I_Upper

    CALL Copy(A%Root,1,MinVec)
    CALL Copy(A%Root,1,MaxVec)

    CALL New(g, 1, MM)
    CALL New(h, 1, MM)
    CALL New(Ax, 1, MM)
    CALL New(Ah, 1, MM)
    CALL New(xOld, 1, MM)
    CALL New(gOld, 1, MM)
    CALL New(hOld, 1, MM)

    CALL spamm_convert_order_2_to_dense (A, A_dense)       

    ALLOCATE (MinVec_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (MaxVec_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (g_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (h_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (Ax_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (Ah_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (xOld_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (gOld_dense(1:SIZE(A_DENSE,1)))
    ALLOCATE (hOld_dense(1:SIZE(A_DENSE,1)))

    call random_number(MinVec_dense)
    call random_number(MaxVec_dense)

    MinVec_dense=MinVec_dense/SQRT(DOT_PRODUCT(MinVec_dense,MinVec_dense))
    MaxVec_dense=MaxVec_dense/SQRT(DOT_PRODUCT(MaxVec_dense,MaxVec_dense))

    DO I=2,2
       IF(I==1)THEN
          x=>MinVec
          x_dense=MinVec_dense
          omega=1d10
       ELSE
          x=>MaxVec
          x_dense=MaxVec_dense
          omega=-1d10
       ENDIF
       xnorm=SpAMM_One/Dot(x,x)
       CALL Multiply(x,xnorm)
       x_dense=x_dense/dot_product(x_dense,x_dense)
      

      ! This call should be redundant  <<<<<<<<<<<>>>>>>>>>>>>>>
      x%Norm=SQRT(Norm(x))
      CALL Multiply(h,SpAMM_Zero)
      CALL Multiply(g,SpAMM_Zero)
      CALL Multiply(xOld,SpAMM_Zero)
      CALL Multiply(hOld,SpAMM_Zero)
      CALL Multiply(gOld,SpAMM_Zero)

      h_dense=0d0
      g_dense=0d0
      xold_dense=0d0
      hold_dense=0d0
      gold_dense=0d0

      
      DO CG=1,NCG
        ! Intermediates
        xx=Dot(x,x)
        CALL Multiply(A%root,x,Ax, Tau)
        xAx=Dot(x,Ax)        
        omega_old=omega
        omega=xAx/xx



        ax_dense=MATMUL(A_dense,x_dense)
        xx_dense=dot_product(x_dense,x_dense)
        xax_dense=dot_product(x_dense,ax_dense)
        omega_dense=xax_dense/dot_product(x_dense,x_dense)

!        WRITE(*,*)xx,xx_dense
!        WRITE(*,*)omega,omega_dense

        IF(I==1)THEN
          RQIMin=omega
          RQIMin_dense=omega_dense
        ELSE
          RQIMax=omega
          RQIMax_dense=omega_dense
        ENDIF
        ! Gradient of extremal quotients (eigenvalues): one is + the other is - ...
        IF(I==1)THEN
          ! g=2*(Ax-omega*x)/xx
          CALL Add(g, SpAMM_Two/xx,Ax,-SpAMM_Two*omega/xx,x)
          g_dense= 2d0*(Ax_dense-omega_dense*x_dense)/xx_dense
        ELSE
          ! g=-2*(Ax-omega*x)/xx
          CALL Add(g,-SpAMM_Two/xx,Ax,SpAMM_Two*omega/xx,x)
          g_dense=-2d0*(Ax_dense-omega_dense*x_dense)/xx_dense
        ENDIF

        WRITE(*,*)'omega = ',omega,omega_old
        IF( SQRT(Dot(g,g)/ABS(Omega)) < CnvrgCrit .AND. CG>16 &
            .OR. (I==1.AND.Omega>Omega_old)                   &
            .OR. (I==2.AND.Omega<Omega_old)                    )EXIT

        IF(CG>1.AND.MOD(CG,15).NE.0)THEN
          !             beta=MAX(SpAMM_Zero,Dot(g,g-gOld)/Dot(gOld,gOld))


          dot_old=Dot(gOld,gOld)
          dot_old_dense=dot_product(gold_dense,gold_dense)

          IF(dot_old.LE.1D-10)THEN
             beta=SpAMM_Zero
             beta_dense=SpAMM_Zero
          ELSE
             beta=MAX(SpAMM_Zero,Dot(g,g)/Dot(gOld,gOld))
             beta_dense=MAX(SpAMM_Zero,dot_product(g_dense,g_dense)/dot_product(gOld_dense,gOld_dense))
          ENDIF
        ELSE
          beta=SpAMM_Zero
          beta_dense=SpAMM_Zero
        ENDIF

        ! h=g+beta*hOld
        CALL Add(h,SpAMM_One,g,beta,hOld)
        h_dense=g_dense+beta_dense*hOld_dense

        h%Norm=SQRT(Norm(h))
        CALL Copy(g,gOld)
        CALL Copy(h,hOld)

        gold_dense=g_dense
        hold_dense=h_dense


        ! Ah=A.h
        CALL Multiply(A%Root,h,Ah, Tau)
        Ah_dense=MATMUL(A_dense,h_dense)

        !-----------------------------------------------------------------------------------------------------------------------
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
        !-----------------------------------------------------------------------------------------------------------------------
        hx_dense =Dot_Product(h_dense,x_dense)
        hh_dense =Dot_Product(h_dense,h_dense)
        xAh_dense=Dot_Product(x_dense,Ah_dense)
        hAx_dense=Dot_Product(h_dense,Ax_dense)
        hAh_dense=Dot_Product(h_dense,Ah_dense)
        ! By symmetry
        xh_dense=hx_dense
        hAx_dense=xAh_dense
        Lambdaplus_Dense=(SpAMM_Two*hh_dense*xAx_dense-SpAMM_Two*hAh_dense*xx_dense+SQRT((-SpAMM_Two*hh_dense*xAx_dense+SpAMM_Two*hAh_dense*xx_dense)**2     &
          -SpAMM_Four*(hAh_dense*hx_dense-SpAMM_Two*hh_dense*xAh_dense+hAh_dense*xh_dense)*(-(hx_dense*xAx_dense)-xAx_dense*xh_dense+SpAMM_Two*xAh_dense*xx_dense)))  &
          /(SpAMM_Two*(hAh_dense*hx_dense-SpAMM_Two*hh_dense*xAh_dense+hAh_dense*xh_dense))
        LambdaMins_Dense=(SpAMM_Two*hh_dense*xAx_dense-SpAMM_Two*hAh_dense*xx_dense-SQRT((-SpAMM_Two*hh_dense*xAx_dense+SpAMM_Two*hAh_dense*xx_dense)**2     &
          -SpAMM_Four*(hAh_dense*hx_dense-SpAMM_Two*hh_dense*xAh_dense+hAh_dense*xh_dense)*(-(hx_dense*xAx_dense)-xAx_dense*xh_dense+SpAMM_Two*xAh_dense*xx_dense)))  &
          /(SpAMM_Two*(hAh_dense*hx_dense-SpAMM_Two*hh_dense*xAh_dense+hAh_dense*xh_dense))
        !
        RQIPlus_dense=(xAx_dense+Lambdaplus_Dense*(xAh_dense+hAx_dense)+hAh_dense*Lambdaplus_Dense**2) &
          /( xx_dense+Lambdaplus_Dense*(xh_dense+hx_dense)  +hh_dense *Lambdaplus_Dense**2)
        RQIMins_dense=(xAx_dense+LambdaMins_Dense*(xAh_dense+hAx_dense)+hAh_dense*LambdaMins_Dense**2) &
          /( xx_dense+LambdaMins_Dense*(xh_dense+hx_dense)  +hh_dense *LambdaMins_Dense**2)

        IF(I==1)THEN
          IF(RQIMins_dense<RQIPlus_dense)THEN
            ! x=x+LambdaMins*h
             x_dense=x_dense+lambdamins_dense*h_dense
          ELSE
            ! x=x+LambdaPlus*h
             x_dense=x_dense+lambdaplus_dense*h_dense
          ENDIF
        ELSE
          IF(RQIMins_dense>RQIPlus_dense)THEN
            ! x=x+LambdaMins*h
             x_dense=x_dense+lambdamins_dense*h_dense
          ELSE
            ! x=x+LambdaPlus*h
             x_dense=x_dense+lambdaplus_dense*h_dense
          ENDIF
        ENDIF

!        x_dense=x_dense/SQRT(dot_product(x_dense,x_dense))
        !-----------------------------------------------------------------------------------------------------------------------

        IF(I==1)THEN
          WRITE(*,33)omega,SQRT(Dot(g,g))/ABS(Omega),CG
!          WRITE(*,33)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG
        ELSE
!          WRITE(*,44)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG
        ENDIF
      END DO

      IF(I==1)THEN
         WRITE(*,33)omega,SQRT(Dot(g,g))/ABS(Omega),CG
         WRITE(*,33)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG
         MinVec=>x
      ELSE
         WRITE(*,44)omega,SQRT(Dot(g,g))/ABS(Omega),CG
         WRITE(*,44)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG         
         MaxVec=>x
      ENDIF

    ENDDO
!    WRITE(*,*)' EIGEN MIn = ',RQIMin,RQIMax
!    WRITE(*,*)' EIGEN MIn = ',RQIMin_dense,RQIMax_dense


33  FORMAT(' MIN E.V. = ',E24.16,', GRAD RQI = ',E16.8,' in ',I4,' NLCG steps')
44  FORMAT(' MAX E.V. = ',E24.16,', GRAD RQI = ',E16.8,' in ',I4,' NLCG steps')

    CALL Delete(g)
    CALL Delete(h)
    CALL Delete(Ax)
    CALL Delete(Ah)
    CALL Delete(xOld)
    CALL Delete(gOld)
    CALL Delete(hOld)

    CALL Delete(MinVec)
    CALL Delete(MaxVec)

    DEALLOCATE(x_dense)
    DEALLOCATE(g_dense)
    DEALLOCATE(h_dense)
    DEALLOCATE(Ax_dense)
    DEALLOCATE(Ah_dense)
    DEALLOCATE(xOld_dense)
    DEALLOCATE(gOld_dense)
    DEALLOCATE(hOld_dense)

  END SUBROUTINE GSOLVE_RQI_Extrema_RQI
END MODULE GSOLVE_RQI

MODULE GSOLVE_INVERSE
  USE SpAMM_algebra
  USE SpAMM_convert
  USE SpAMM_management
  USE SpAMM_types
  USE SpAMM_utilities
  USE GSOLVE_RQI
  IMPLICIT NONE

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1

CONTAINS

  FUNCTION GSOLVE_SPECTRAL_Shift( X, low_prev, high_prev, low_new, high_new ) 
    TYPE(spamm_matrix_order_2) , POINTER                   :: GSOLVE_SPECTRAL_Shift
    TYPE(spamm_matrix_order_2) , POINTER, INTENT(INOUT)    :: X
    REAL(SpAMM_KIND), OPTIONAL,INTENT(IN)    :: low_prev, high_prev, low_new, high_new   
    INTEGER                                  :: M,N
    REAL(SpAMM_KIND)                         :: SHFT,SCAL 

    M=X%M
    N=X%N
    !!!!!    shft=low_new+(x-low_prev)*(high_new-low_new)/(high_prev-low_prev)
    SHFT=low_new-low_prev*(high_new-low_new)/(high_prev-low_prev)
    SCAL=(high_new-low_new)/(high_prev-low_prev)

!    WRITE(*,*)' prev low = ',low_prev,' new  low = ',low_new

    CALL Multiply( X, SCAL )
    CALL Add(X%Root, SHFT, M, N)
    GSOLVE_SPECTRAL_Shift => X ! dear me, is that kosher?
  END FUNCTION GSOLVE_SPECTRAL_Shift

  FUNCTION GSOLVE_SPECTRAL_InvSqrt_Scaled_NS( X, sc ) 
    TYPE(spamm_matrix_order_2) , POINTER                   :: GSOLVE_SPECTRAL_InvSqrt_Scaled_NS
    TYPE(spamm_matrix_order_2) , POINTER, INTENT(INOUT)    :: X
    REAL(SpAMM_KIND),      INTENT(IN)      :: sc 
    REAL(SpAMM_KIND)                       :: SHFT,SCAL 
    INTEGER                                :: M,N
    M=X%M
    N=X%N
    SHFT=5d-1*SQRT(sc)*(3.0d0)
    SCAL=5d-1*SQRT(sc)*(  -sc)
    CALL Multiply( X, SCAL )
    CALL Add(X%Root, SHFT, M, N)
    GSOLVE_SPECTRAL_InvSqrt_Scaled_NS => X ! dear me, is that kosher?
  END FUNCTION GSOLVE_SPECTRAL_InvSqrt_Scaled_NS

  REAL(SpAMM_KIND) FUNCTION ScaleInvSqrt(xo)
    REAL(SpAMM_KIND) :: xo
    ScaleInvSqrt=MIN( Approx3, 3.d0/(1d0+SQRT(xo)+xo) )    
  END FUNCTION ScaleInvSqrt

  SUBROUTINE GSOLVE_NORMALIZE_MATRIX(S, EvMax )
    TYPE(spamm_matrix_order_2) , POINTER, INTENT(INOUT) :: S
    TYPE(BiTree) , POINTER                :: MinV=>NULL(),MaxV=>NULL()
    INTEGER                               :: i,j,M,N,M_p,N_p
    REAL(SpAMM_KIND)                      :: Tau = 1d-6
    REAL(SpAMM_KIND)                      :: EvMin, EvMax, sc, xo, xo_prev, xo_analytic, xmax,delta, FillN

    real(spamm_kind), dimension(:, :), allocatable :: X_dense,Z_Dense,S_Dense

    integer :: LWORK
    integer :: LIWORK
    real(kind(0d0)), allocatable :: eval(:), SHalf(:,:),SHlfI(:,:),evec(:, :), work(:)
    integer, allocatable :: iwork(:)
    integer :: info

    interface
       subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
         CHARACTER          JOBZ, UPLO
         INTEGER            INFO, LDA, LIWORK, LWORK, N
         INTEGER            IWORK( * )
         DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
       end subroutine dsyevd
    end interface


    M=S%M
    N=S%N
    M_p=S%Root%I_upper-S%Root%I_lower+1
    N_p=S%Root%J_upper-S%Root%J_lower+1


    CALL New(MinV, 1, N_p)
    CALL New(MaxV, 1, N_p)
    CALL GSOLVE_RQI_Extrema_RQI(S,EvMin,EvMax,MinV,MaxV,Tau,Tau)
    WRITE(*,*)' MinV = ',EvMin,', MaxV = ',EvMax

    CALL Multiply(S,1d0/EvMax)

  END SUBROUTINE GSOLVE_NORMALIZE_MATRIX


  SUBROUTINE GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(S, Z, Tau)

    TYPE(spamm_matrix_order_2) , POINTER, INTENT(IN)    :: S
    TYPE(spamm_matrix_order_2) , POINTER, INTENT(INOUT) :: Z
    TYPE(spamm_matrix_order_2) , POINTER                :: X =>NULL(),T1=>NULL(),T2=>NULL(), &
                                                           X2=>NULL(),  ZT=>NULL()

    INTEGER                               :: i,j,M,N,M_p,N_p
    REAL(SpAMM_KIND),       INTENT(IN)    :: Tau
    REAL(SpAMM_KIND)                      :: EvMin, EvMax, sc, xo, xo_prev, xo_analytic, xmax,delta, FillN, FillN_prev

    real(spamm_kind), dimension(:, :), allocatable :: X_dense,Z_Dense,S_Dense

    integer :: LWORK
    integer :: LIWORK
    real(kind(0d0)), allocatable :: eval(:), SHalf(:,:),SHlfI(:,:),evec(:, :), work(:)
    integer, allocatable :: iwork(:)
    integer :: info

    REAL(SpAMM_KIND) :: EvMinDense,EvMinSpAMM,EvMinShift,EvMinMapped, &
                       EvMaxDense,EvMaxSpAMM,EvMaxShift,EvMaxMapped, XDiffMax

    interface
       subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
         CHARACTER          JOBZ, UPLO
         INTEGER            INFO, LDA, LIWORK, LWORK, N
         INTEGER            IWORK( * )
         DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
       end subroutine dsyevd
    end interface

    M=S%M
    N=S%N
    M_p=S%Root%I_upper-S%Root%I_lower+1
    N_p=S%Root%J_upper-S%Root%J_lower+1
    !------------------------------------ --------------------------------------------
    CALL spamm_convert_order_2_to_dense (S, X_dense)       
    CALL spamm_convert_order_2_to_dense (S, S_dense)       
    N = SIZE(X_dense,1)
    LWORK = 1+6*N+2*N**2
    LIWORK = 3+5*N    
    allocate(eval(N))
    allocate(work(LWORK))
    allocate(iwork(LIWORK))
    call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
    CALL spamm_convert_order_2_to_dense (S, X_dense)       
    xo=eval(1)
    xo_analytic=xo
    !------------------------------------ --------------------------------------------
    Z => SpAMM_identity_matrix(M,N)
    ZT => SpAMM_identity_matrix(M,N)
    X => SpAMM_identity_matrix(M,N)
    X2 => SpAMM_identity_matrix(M,N)
    T1 => SpAMM_identity_matrix(M,N)
    T2 => SpAMM_identity_matrix(M,N)
    !------------------------------------ --------------------------------------------
    !
    WRITE(*,22)
22  FORMAT(' It,      Tau,  ~Zt.S.Z, MinExact, MinSpAMM, MinShift, MinMappd,', &
           ' MaxExact, MaxSpAMM, MaxShift, MaxMappd,  Idealty, FillIn')
    !
    FillN=1d10
    !
    DO i = 1, 22
       !
       CALL spamm_convert_order_2_to_dense (Z, Z_dense)       
       ZT=>SpAMM_convert_dense_to_matrix_2nd_order(TRANSPOSE(Z_dense))
       !       CALL SpAMM_Transpose_QuTree( Z%root, ZT%root )

       CALL Multiply( T1, SpAMM_Zero)       
       CALL Multiply( ZT, S,  T1, 1d-3*Tau )
       CALL Multiply( T1, Z, X, Tau )
       !--------------------------------------------------------------------------------
       ! Debug ....
       X_DENSE=MATMUL(TRANSPOSE(Z_DENSE),MATMUL(S_DENSE,Z_DENSE))
       X2=>SpAMM_convert_dense_to_matrix_2nd_order(X_dense)
       CALL Add(X2, X, -1d0, 1d0)
       XDiffMax=absmax(X2)
       call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
       EvMinDense=eval(1)
       EvMaxDense=eval(N)
       CALL spamm_convert_order_2_to_dense (X, X_dense)       
       call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
       EvMinSpAMM=eval(1)
       EvMaxSpAMM=eval(N)
       ! end debug
       !--------------------------------------------------------------------------------
       ! No, we don't compute the lowest eigenvalue in practice.  Deal w/it...
       sc=ScaleInvSqrt(0d0)
       ! instead, lets monitor the trace for convergence:
       FillN_prev=FillN
       FillN=DBLE(N-Trace(X))/DBLE(N)       
       !        
       IF(FillN>0.4d0)THEN
          delta=1d-1
          X => GSOLVE_SPECTRAL_Shift( X, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
       ELSE
          sc=1d0
       ENDIF
       !--------------------------------------------------------------------------------
       ! Debug ....
       CALL spamm_convert_order_2_to_dense (X, X_dense)       
       call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
       EvMinShift=eval(1)
       EvMaxShift=eval(N)
       ! end debug
       !--------------------------------------------------------------------------------
       X =>  GSOLVE_SPECTRAL_InvSqrt_Scaled_NS( X, sc ) 

       CALL Multiply(T2, SpAMM_Zero)
       CALL multiply( Z, X, T2, Tau ) 
       CALL Copy(T2,Z)
       !--------------------------------------------------------------------------------
       ! Debug ....
       CALL spamm_convert_order_2_to_dense (X, X_dense)       
       call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
       EvMinMapped=eval(1)
       EvMaxMapped=eval(N)

       WRITE(*,33)I,Tau,XDiffMax,EvMinDense,EvMinSpAMM,EvMinShift,EvMinMapped, &
                  EvMaxDense,EvMaxSpAMM,EvMaxShift,EvMaxMapped, &
                  EvMinDense/xo_analytic, FillN

33     FORMAT(I3,', ',12(E8.2,', '))
       ! end debug
       !--------------------------------------------------------------------------------

       IF(FillN<0.1d0.AND.FillN>FillN_prev)THEN
          RETURN
       ELSEIF(FillN<Tau)THEN
          RETURN
       ENDIF
       ! best acceleration we can hope for
       xo_analytic=xo_analytic*(9d0/4d0)*sc
       !
    ENDDO

    call delete(X)
    call delete(X2)
    call delete(ZT)
    call delete(T1)
    call delete(T2)

    deallocate(eval)
    deallocate(work)
    deallocate(iwork)

    deallocate(X_Dense)
    deallocate(Z_Dense)
    deallocate(S_Dense)

  END SUBROUTINE  GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT


  SUBROUTINE SetScalarSpectrum(CondS,S,N)
    INTEGER :: M,I,N
    REAL*8 :: CondS
    REAL*8, DIMENSION(N) :: S
    integer,parameter :: seed = 86456
    REAL*8 :: low1,low2,high1,high2, dx, XMin,XMax

    call srand(seed)
    !    
    low1 = 1d10
    high1=-1d10
    DO i=1,N
       s(i)=RAND()
       low1 =MIN(low1, s(i))
       high1=MAX(high1,s(i))
    ENDDO         
    !
    low2=-LOG10(CondS)    
    high2=0

    DO i=1,N
       s(i)=shft(s(i), dx, low_old=low1, high_old=high1, low_new=low2, high_new=high2)
       s(i)=10d0**s(i)
    ENDDO

  END SUBROUTINE SetScalarSpectrum

  REAL*8 function shft(x, dx, low_old, high_old, low_new, high_new)
    REAL*8, INTENT(INOUT)  :: x, dx
    REAL*8, OPTIONAL :: low_old, high_old, low_new, high_new    
    dx=low_new-low_old
    shft=low_new+(x-low_old)*(high_new-low_new)/(high_old-low_old)
  END function shft

END MODULE GSOLVE_INVERSE

program SpAMM_lowdin

  USE SpAMMpack
  USE SpAMM_PROJECT
  USE test_utilities

  USE GSOLVE_INVERSE

  implicit none

#ifdef LAPACK_FOUND
  interface
     subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LIWORK, LWORK, N
       INTEGER            IWORK( * )
       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
     end subroutine dsyevd
  end interface
#endif

  ! COND(5)
  real(kind(0d0)), parameter ::   Tau_1=1d-4
  real(kind(0d0)), parameter ::   Tau_2=1d-8
  real(kind(0d0)), parameter ::   Tau_3=1d-12

  type(SpAMM_matrix_order_2), pointer :: S => null(), X => null(), X1=>NULL(), X2=>NULL(), &
       ZT=>NULL(),Z => null(), Z1=>NULL(),Z2 => null(), Id => null(),T1=>NULL()

  real(kind(0d0)), allocatable :: x_dense(:,:),S_dense(:, :), Id_dense(:, :), SNew(:,:),Z_dense(:,:)
  character(len = 1000) :: matrix_filename
  real :: start_time, end_time
  integer :: i
  real(kind(0d0)) :: max_diff, evmax

  integer :: M,N
  integer :: LWORK
  integer :: LIWORK
  real(kind(0d0)), allocatable :: eval(:), SHalf(:,:),SHlfI(:,:),evec(:, :), work(:)
  integer, allocatable :: iwork(:)
  integer :: info
  logical :: fakeit


  if(command_argument_count() == 0) then
     write(*, "(A)") "Please specify the matrix file to USE (the overlap matrix)"
     error stop
  else if(command_argument_count() > 1) then
     write(*, "(A)") "Too many command line arguments. One is plenty."
     error stop
  else
     call get_command_argument(1, matrix_filename)
  endif

  call read_MM(matrix_filename, S_dense)

  N = size(S_dense, 1)
  LWORK = 1+6*N+2*N**2
  LIWORK = 3+5*N

  allocate(eval(N))
  allocate(work(LWORK))
  allocate(iwork(LIWORK))
  allocate(snew(1:N,1:N))

!  FakeIt=.TRUE.
  FakeIt=.FALSE.
  IF(FakeIt)THEN
     evec=S_dense
     call dsyevd("V", "U", N, evec, N, eval, work, LWORK, iwork, LIWORK, info)
     CALL SetScalarSpectrum(1d10,eval,N)
     DO i=1,N
        S_dense(:,I)=evec(:,i)*eval(i)
     ENDDO
     SNew=MATMUL(S_dense,TRANSPOSE(evec))
     S => SpAMM_convert_dense_to_matrix_2nd_order(SNew)
  ELSE
     S => SpAMM_convert_dense_to_matrix_2nd_order(S_dense)
     CALL GSOLVE_NORMALIZE_MATRIX(S, EvMax )
  ENDIF

  CALL Copy(S,X1)
  CALL GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(X1, Z1, Tau_1 )

  M=S%M
  N=S%N
  ZT => SpAMM_identity_matrix(M,N)
  T1 => SpAMM_identity_matrix(M,N)
  CALL Multiply( ZT, SpAMM_Zero)       
  CALL Multiply( T1, SpAMM_Zero)       

  CALL SpAMM_Transpose_QuTree( Z1%root, ZT%root )

  CALL Copy(S,X2)

  CALL spamm_convert_order_2_to_dense (Z1, Z_dense)       
  ZT=>SpAMM_convert_dense_to_matrix_2nd_order(TRANSPOSE(Z_dense))

  CALL Multiply( X2,  Z1, T1, Tau_2 )
  CALL Multiply( ZT, T1,  X2, Tau_2      )

  CALL GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(X2, Z2, Tau_2 )

  STOP

111 CONTINUE


  evec=S_dense
  CALL SetScalarSpectrum(1d10,eval,N)
  DO i=1,N
     S_dense(:,I)=S_dense(:,i)*eval(i)
  ENDDO
  SNew=MATMUL(evec,TRANSPOSE(S_dense))
  S => SpAMM_convert_dense_to_matrix_2nd_order(SNew)

  call cpu_time(start_time)
  call GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(S, Z, 1d-3)!tolerance) 
  call cpu_time(end_time)

  STOP  

end program SpAMM_lowdin
