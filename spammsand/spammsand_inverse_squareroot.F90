#define DENSE_DIAGNOSTICS
! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan  -Wextra -std=f2008 "

module SpAMMsand_inverse_squareroot
  USE spammpack 

  implicit none

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1

#ifdef DENSE_DIAGNOSTICS
  real(spamm_kind), dimension(:,:), ALLOCATABLE ::   i_d, s_d,                                & 
     z_k_stab,     z_k1_stab,     x_k_stab,     x_k1_stab,      m_x_k1_stab,                   &   
     z_tld_k_stab, z_tld_k1_stab, x_tld_k_stab, x_tld_k1_stab,  m_x_tld_k1_stab,               &
     z_k_dual,     z_k1_dual,     x_k_dual,     x_k1_dual,      m_x_k1_dual,      y_k_dual,     y_k1_dual,      &
     z_tld_k_dual, z_tld_k1_dual, x_tld_k_dual, x_tld_k1_dual,  m_x_tld_k1_dual,  y_tld_k_dual, y_tld_k1_dual  

  REAL(spamm_kind) :: scal_shift_dual,shft_shift_dual, scal_mapp_dual, shft_mapp_dual
  REAL(spamm_kind) :: scal_shift_stab,shft_shift_stab, scal_mapp_stab, shft_mapp_stab
  
  integer :: LWORK
  integer :: LIWORK
  real(kind(0d0)), allocatable :: eval(:), SHalf(:,:),SHlfI(:,:),evec(:, :), work(:)
  integer, allocatable :: iwork(:)
  integer :: info, N  

  interface
     subroutine dsyevd ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LIWORK, LWORK, N
       INTEGER            IWORK( * )
       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
     end subroutine dsyevd
  end interface

#endif
  
contains

  FUNCTION IntToChar(I)
    CHARACTER(LEN=*), PARAMETER :: INTERNAL_INT_FMT='(I22)'
    INTEGER  :: I
    CHARACTER(LEN=22) :: IntToChar
    WRITE(UNIT=IntToChar,FMT=INTERNAL_INT_FMT)I
    IntToChar=ADJUSTL(IntToChar)
  END FUNCTION IntToChar


  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, tau, DoDuals, second, kount)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN) :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER             :: x,z
    TYPE(SpAMM_tree_2d_symm) , POINTER             :: x_stab, z_stab, x_dual, z_dual, y_dual, y_tmp, z_tmp

    REAL(SpAMM_KIND) :: Tau
!    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau
!    LOGICAL, INTENT(IN)                              :: first
    LOGICAL                                           :: DoDuals,second
    INTEGER, INTENT(INOUT)                            :: kount
    INTEGER                                           :: i,  j, k
    REAL(SpAMM_KIND)                                  :: sc, TrX, tau_xtra
    REAL(SpAMM_KIND)                                  :: x_work,z_work,y_work,zs_work,sz_work
    REAL(SpAMM_KIND)                                  :: x_fill,z_fill,y_fill,zs_fill,sz_fill,y_norm,zs_norm
    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev
    REAL(SpAMM_KIND)                                  :: tau_zdotz,sz_norm



    

!    tau_xtra=1d-2*tau

#ifdef DENSE_DIAGNOSTICS
#else
    !
    IF(DoDuals)tHeN
#endif
       ! y_0 => s
       y_dual => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       y_dual => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_o = y_dual, threshold_O = SpAMM_normclean ) 
       
       !  z_0 => I
       z_dual => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       z_dual => SpAMM_set_identity_2d_symm( s%frill%ndimn, in_o = z_dual )

       !  x_0 => s
       x_dual => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       x_dual => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_O = x_dual, threshold_O = SpAMM_normclean ) 

       WRITE(*,*)'d tr = ',SpAMM_trace_tree_2d_symm_recur(x_dual)

#ifdef DENSE_DIAGNOSTICS
#else
    ELSE
#endif
       !  z_0 => I
       z_stab => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       z_stab => SpAMM_set_identity_2d_symm( s%frill%ndimn, in_o = z_stab )

       !  x_0 => s
       x_stab => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       x_stab => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_O = x_stab, threshold_O = SpAMM_normclean ) 

       WRITE(*,*)'s tr = ',SpAMM_trace_tree_2d_symm_recur(x_stab)

#ifdef DENSE_DIAGNOSTICS
#else
    endIF
#endif

    y_tmp  => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
    z_tmp  => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )

#ifdef DENSE_DIAGNOSTICS

    n=s%frill%ndimn(1)

    ALLOCATE(i_d     (1:N,1:N))
    ALLOCATE(s_d     (1:N,1:N))

    allocate(z_k_stab     (1:n,1:n));  allocate(z_k1_stab     (1:n,1:n));  allocate(x_k_stab    (1:n,1:n)); 
    allocate(x_k1_stab    (1:n,1:n));  allocate(m_x_k1_stab   (1:n,1:n))        

    allocate(z_tld_k_stab (1:n,1:n)); allocate(z_tld_k1_stab  (1:n,1:n));  allocate(x_tld_k_stab(1:n,1:n)); 
    allocate(x_tld_k1_stab(1:n,1:n)); allocate(m_x_tld_k1_stab(1:n,1:n))

    allocate(z_k_dual     (1:n,1:n)); allocate(z_k1_dual      (1:n,1:n));  allocate(x_k_dual    (1:n,1:n)); 
    allocate(x_k1_dual    (1:n,1:n)); allocate(m_x_k1_dual    (1:n,1:n));  allocate(y_k_dual    (1:n,1:n));  allocate(y_k1_dual    (1:n,1:n))       

    allocate(z_tld_k_dual (1:n,1:n)); allocate(z_tld_k1_dual  (1:n,1:n));  allocate(x_tld_k_dual(1:n,1:n));  
    allocate(x_tld_k1_dual(1:n,1:n)); allocate(m_x_tld_k1_dual (1:n,1:n)); allocate(y_tld_k_dual(1:n,1:n));  allocate(y_tld_k1_dual(1:n,1:n)) 

    ! I=diag(1)
    i_d=SpAMM_zero
    DO k=1,n
       i_d(k,k)=SpAMM_one
    ENDDO
    !
    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d  )
    ! --------------------
    z_k1_stab    =i_d
    z_tld_k1_stab=i_d
    x_k1_stab    =s_d
    x_tld_k1_stab=s_d
    ! --------------------
    x_k1_dual    =s_d
    x_tld_k1_dual=s_d
    z_k1_dual    =i_d
    z_tld_k1_dual=i_d
    y_k1_dual    =s_d
    y_tld_k1_dual=s_d
    ! --------------------
    DoDuals=.TRUE.

#endif
    kount=1
    sz_norm=1d0
    FillN=1d10
    DO i = 0, 17

       tau=1d-10
       tau_xtra=1d-10

#ifdef DENSE_DIAGNOSTICS
       WRITE(*,*)' n = ',I, ' tr = ',SpAMM_trace_tree_2d_symm_recur(x_dual),SpAMM_trace_tree_2d_symm_recur(x_stab)
#else
       ! check the trace for convergence:
       FillN_prev=FillN
       IF(DoDuals)THEN
          TrX=SpAMM_trace_tree_2d_symm_recur(x_dual)
       ELSE
          TrX=SpAMM_trace_tree_2d_symm_recur(x_stab)
       ENDIF
       FillN = abs( dble(s%frill%ndimn(1)) - TrX )/dble(s%frill%ndimn(1))       
       IF(DoDuals)THEN
          WRITE(*,33)tau, kount, TrX, FillN, y_work*1d2 , z_work*1d2 , x_work*1d2 
       elsE
          WRITE(*,34)tau, kount, TrX, FillN, y_work*1d2 , z_work*1d2 , x_work*1d2 
       ENDIF
#endif

33     format('  dual Tr< ',e6.1,', n=',i2,' > = ', F18.10,' dN=',e10.3,', y_wrk: ',f10.5,'%, z_wrk:',f10.5,'%, x_wrk:',f10.5,'%')
34     format('  stab Tr< ',e6.1,', n=',i2,' > = ', F18.10,' dN=',e10.3,', y_wrk: ',f10.5,'%, z_wrk:',f10.5,'%, x_wrk:',f10.5,'%')

       IF(i>2 .and. FillN<0.1d0 .AND. FillN>FillN_prev )then

!          WRITE(*,*)' fill n = ',filln,' filln_prev ',filln_prev
!          write(*,*)' elevation' 

!          RETURN  ! Elevation

       endif

       IF( FillN <  Tau )then

!          WRITE(*,*)' fill n = ',filln,' tau**2 = ',tau
!          write(*,*)' anhiliaiton '

!          RETURN  ! Anihilation

       end IF

#ifdef DENSE_DIAGNOSTICS

       m_x_k1_stab=x_k1_stab
       m_x_k1_dual=x_k1_dual

       m_x_tld_k1_stab=x_tld_k1_stab
       m_x_tld_k1_dual=x_tld_k1_dual

       scal_shift_stab=1d0
       shft_shift_stab=0d0 
       scal_mapp_stab =1d0
       shft_mapp_stab =0d0

       scal_shift_dual=1d0
       shft_shift_dual=0d0 
       scal_mapp_dual =1d0
       shft_mapp_dual =0d0

#endif

!!$       IF(DoDuals)THEN
!!$          sc=1d0
!!$
!!$          write(*,*)' a '
!!$       ELSE
          
          IF(FillN>0.4d0)THEN

             delta=8.d-2 ! maybe this should be a variable too, passed in?

!!$#ifdef DENSE_DIAGNOSTICS
!!$             x_stab => spammsand_shift_tree_2d( x_stab, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta , stab= .TRUE. )
!!$             x_dual => spammsand_shift_tree_2d( x_dual, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta , stab= .FALSE. )
!!$#else
!!$             IF(DoDuals)THEN
!!$                x_dual => spammsand_shift_tree_2d( x_dual, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
!!$             ELSE
!!$                x_stab => spammsand_shift_tree_2d( x_stab, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
!!$             ENDIF
!!$#endif
!!$             sc=spammsand_scaling_invsqrt(SpAMM_zero)
!!$           ELSEIF(FillN>0.1d0)THEN
!!$             delta=1.0d-2 ! maybe this should be a variable too, passed in?
!!$             x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
!!$             sc=spammsand_scaling_invsqrt(SpAMM_half)

          ELSE

             sc=1d0

          ENDIF
 
!       ENDIF

       sc=1d0

       WRITE(*,*)' sc = ',sc
       
#ifdef DENSE_DIAGNOSTICS       

       x_stab => spammsand_scaled_invsqrt_mapping( x_stab, sc, .TRUE.  )
       z_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_stab, x_stab, tau , nt_O=.TRUE., in_O = z_tmp , stream_file_O='z_stab'//inttoCHAR(I) )
       z_stab => SpAMM_tree_2d_symm_copy_tree_2d_symm( z_tmp, in_O = z_stab, threshold_O = tau ) 

       x_dual => spammsand_scaled_invsqrt_mapping( x_dual, sc, .FALSE. )
       z_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_dual, x_dual, tau , nt_O=.TRUE., in_O = z_tmp , stream_file_O='z_dual'//inttoCHAR(I) )
       z_dual=> SpAMM_tree_2d_symm_copy_tree_2d_symm( z_tmp, in_O = z_dual, threshold_O = tau ) 

#else
       IF(DoDuals)THEN
          ! m[x_n-1,c]   
          x_dual => spammsand_scaled_invsqrt_mapping( x_dual, sc)
          ! |z_n> =  <z_n-1| m[x_n-1]  
          z_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_dual, x_dual, tau , nt_O=.TRUE., in_O = z_tmp , stream_file_O='z_dual'//inttoCHAR(I) )
          z_dual=> SpAMM_tree_2d_symm_copy_tree_2d_symm( z_tmp, in_O = z_dual, threshold_O = tau ) 
          z_work=z_tmp%frill%flops/dble(z_tmp%frill%ndimn(1))**3
          z_fill=z%frill%non0s
       ELSE 
          ! m[x_n-1,c]   
          x_stab => spammsand_scaled_invsqrt_mapping( x_stab, sc)
          ! |z_n> =  <z_n-1| m[x_n-1]  
          z_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_stab, x_stab, tau , nt_O=.TRUE., in_O = z_tmp , stream_file_O='z_stab'//inttoCHAR(I) )
          z_stab=> SpAMM_tree_2d_symm_copy_tree_2d_symm(  z_tmp, in_O = z_stab, threshold_O = tau ) 
          z_work=z_tmp%frill%flops/dble(z_tmp%frill%ndimn(1))**3
          z_fill=z%frill%non0s
       ENDIF
#endif
 
#ifdef DENSE_DIAGNOSTICS
       CALL SpAMM_convert_tree_2d_symm_to_dense( z_stab, z_tld_k_stab )
       CALL SpAMM_convert_tree_2d_symm_to_dense( z_dual, z_tld_k_dual )

       z_k_stab = MATMUL( z_k1_stab, m_x_k1_stab )
       z_k_dual = MATMUL( z_k1_dual, m_x_k1_dual )

#else
       IF(DoDuals)tHeN
#endif
          ! |y_n> = m[x_n-1]|y_n-1>
          y_tmp  => SpAMM_tree_2d_symm_times_tree_2d_symm( x_dual, y_dual, tau_xtra , nt_O=.TRUE., in_O = y_tmp , stream_file_O='y_dual_'//inttoCHAR(I) )
          y_dual => SpAMM_tree_2d_symm_copy_tree_2d_symm( y_tmp, in_O = y_dual, threshold_O = tau_xtra )

          y_work=y_tmp%frill%flops/dble(y_tmp%frill%ndimn(1))**3
          y_norm=SQRT(y_tmp%frill%norm2)

          ! |x_n> = <y_n|z_n>
          x_dual => SpAMM_tree_2d_symm_times_tree_2d_symm( y_dual, z_dual, tau   , nt_O=.TRUE. , in_O = x_dual )
          x_work=x_dual%frill%flops/dble(x_dual%frill%ndimn(1))**3

#ifdef DENSE_DIAGNOSTICS
          y_k_dual = MATMUL( m_x_k1_dual,  y_k1_dual )
          CALL SpAMM_convert_tree_2d_symm_to_dense( y_dual, y_tld_k_dual )
          x_k_dual = MATMUL( y_k_dual , z_k_dual )
          CALL SpAMM_convert_tree_2d_symm_to_dense( x_dual, x_tld_k_dual )
#else
       ELSE
#endif
          ! | y_n> = <zt_n|s>
          y_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_stab, s,  tau_XTRA  , NT_O=.FALSE. , in_O = y_tmp )
          y_work=y_tmp%frill%flops/dble(y_tmp%frill%ndimn(1))**3
          y_norm=SQRT(y_tmp%frill%norm2)
          
          ! |x_n> = <y_n|z_n>
          x_stab => SpAMM_tree_2d_symm_times_tree_2d_symm( y_tmp, z_stab,  tau , NT_O=.TRUE. , in_O = x_stab , stream_file_O='x_stab_'//inttoCHAR(I) )
          x_work=x_stab%frill%flops/dble(x_stab%frill%ndimn(1))**3
          
#ifdef DENSE_DIAGNOSTICS

          CALL SpAMM_convert_tree_2d_symm_to_dense( x_stab, x_tld_k_stab )          
          x_k_stab = MATMUL( TRANSPOSE( z_k_stab ) , MATMUL( s_d, z_k_stab ) )

          CALL SpAMMsand_Error_Analysis(kount, tau)

          y_k1_dual=y_k_dual
          z_k1_dual=z_k_dual
          x_k1_dual=x_k_dual
          z_k1_stab=z_k_stab
          x_k1_stab=x_k_stab

          y_k1_dual=y_k_dual
          z_k1_dual=z_k_dual
          x_k1_dual=x_k_dual
          z_k1_stab=z_k_stab
          x_k1_stab=x_k_stab

          y_tld_k1_dual=y_tld_k_dual
          z_tld_k1_dual=z_tld_k_dual
          x_tld_k1_dual=x_tld_k_dual
          z_tld_k1_stab=z_tld_k_stab
          x_tld_k1_stab=x_tld_k_stab
#else
       endif
#endif
       kount=kount+1
    END DO
#ifdef DENSE_DIAGNOSTICS
    STOP ' This has been a diagnostics run, please see fort.99 '
#endif

!    CaLl SpAMM_destruct_tree_2d_symm_recuR (Y)
!    CaLl SpAMM_destruct_tree_2d_symm_recuR (Y)
!    CaLl SpAMM_destruct_tree_2d_symm_recuR (Y_tmp)
!    CaLl SpAMM_destruct_tree_2d_symm_recuR (Z_tmp)
   
  END SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot

#ifdef DENSE_DIAGNOSTICS

  SUBROUTINE SpAMMsand_Error_Analysis(kount, tau)

    real(spamm_kind) :: Tau, FillN_stab, FillN_dual

    integer :: i,kount,j,k

    real(spamm_kind), dimension(:,:), ALLOCATABLE :: &
    dz_stab, dz_hat_stab, dx_stab, dx_hat_stab,       &
    dy_dual, dy_hat_dual,dz_dual, dz_hat_dual,       &
    dx_dual, dx_hat_dual,mp_stab, mp_dual,           &
    dz_k_dx_k1_stab,dz_k_dz_k1_stab,dz_k_dz_k1_dual, &
    dy_k_dy_k1_dual,dz_k_dx_k1_dual,dy_k_dx_k1_dual, &
    x_tld_k_of_xk1_stab,x_tld_k_of_zk1_stab,         &
    x_tld_k_of_xk1_dual,x_tld_k_of_yk1_dual,x_tld_k_of_zk1_dual,       &
    xp_tld_k_stab_gateaux,zp_tld_k_stab_gateaux,xp_tld_k_dual_gateaux, & 
    yp_tld_k_dual_gateaux,zp_tld_k_dual_gateaux,                       &
    FpX_stab,FpZ_stab,FpX_dual,FpZ_dual,FpY_dual  

    allocate(dz_stab                (1:N,1:N) )
    allocate(dz_hat_stab            (1:N,1:N) )
    allocate(dx_stab                (1:N,1:N) )
    allocate(dx_hat_stab            (1:N,1:N) )

    allocate(dy_dual                (1:N,1:N) )
    allocate(dy_hat_dual            (1:N,1:N) )
    allocate(dz_dual                (1:N,1:N) )
    allocate(dz_hat_dual            (1:N,1:N) )
    allocate(dx_dual                (1:N,1:N) )
    allocate(dx_hat_dual            (1:N,1:N) )
    allocate(mp_stab                (1:N,1:N) )
    allocate(mp_dual                (1:N,1:N) )

    allocate(dz_k_dx_k1_stab        (1:N,1:N) )
    allocate(dz_k_dz_k1_stab        (1:N,1:N) )

    allocate(dz_k_dz_k1_dual        (1:N,1:N) )
    allocate(dy_k_dy_k1_dual        (1:N,1:N) )

    allocate(dz_k_dx_k1_dual        (1:N,1:N) )
    allocate(dy_k_dx_k1_dual        (1:N,1:N) )

    allocate(x_tld_k_of_xk1_stab    (1:N,1:N) )
    allocate(x_tld_k_of_zk1_stab    (1:N,1:N) )
    allocate(x_tld_k_of_xk1_dual    (1:N,1:N) )
    allocate(x_tld_k_of_yk1_dual    (1:N,1:N) )
    allocate(x_tld_k_of_zk1_dual    (1:N,1:N) )
    allocate(xp_tld_k_stab_gateaux  (1:N,1:N) )
    allocate(zp_tld_k_stab_gateaux  (1:N,1:N) )
    allocate(xp_tld_k_dual_gateaux  (1:N,1:N) )
    allocate(yp_tld_k_dual_gateaux  (1:N,1:N) )
    allocate(zp_tld_k_dual_gateaux  (1:N,1:N) )
    allocate(FpX_stab               (1:N,1:N) )
    allocate(FpZ_stab               (1:N,1:N) )
    allocate(FpX_dual               (1:N,1:N) )
    allocate(FpZ_dual               (1:N,1:N) )
    allocate(FpY_dual               (1:N,1:N) )

    !-----------------

    WRITE(*,*)SUM(z_tld_k1_stab**2) ,SUM(z_k1_stab**2)
 
    dz_stab     = z_tld_k1_stab - z_k1_stab
    dz_hat_stab = dz_stab/SQRT(SUM(dz_stab**2))
    WRITE(*,*)' dz_stab ',SQRT(SUM(dz_stab**2))

    dx_stab     = x_tld_k1_stab - x_k1_stab
    dx_hat_stab = dx_stab/SQRT(SUM(dx_stab**2))
    WRITE(*,*)' dx_stab ',SQRT(SUM(dx_stab**2))

    dy_dual     = y_tld_k1_dual - y_k1_dual
    dy_hat_dual = dy_dual/SQRT(SUM(dy_dual**2))
    WRITE(*,*)' dy_dual ',SQRT(SUM(dy_dual**2))

    dz_dual     = z_tld_k1_dual - z_k1_dual
    dz_hat_dual = dz_dual/SQRT(SUM(dz_dual**2))
    WRITE(*,*)' dz_dual ',SQRT(SUM(dz_dual**2))


    dx_dual     = x_tld_k1_dual - x_k1_dual
    dx_hat_dual = dx_dual/SQRT(SUM(dx_dual**2))

    mp_stab     = scal_mapp_stab * scal_shift_stab * dx_hat_stab
    mp_dual     = scal_mapp_dual * scal_shift_dual * dx_hat_dual

    !-----------------

    dz_k_dx_k1_stab = MATMUL( z_tld_k1_stab, mp_stab ) 
    dz_k_dz_k1_stab = MATMUL( dz_hat_stab, m_x_tld_k1_stab )



    WRITE(*,*) ' dz_hat_dual, m_x_tld_k1_dual ',SUM(dz_hat_dual**2),SUM(m_x_tld_k1_dual**2)

    dy_k_dy_k1_dual = MATMUL( m_x_tld_k1_dual, dy_hat_dual )




    x_tld_k_of_xk1_stab   = MATMUL( MATMUL( TRANSPOSE(MATMUL(z_tld_k1_stab , m_x_k1_stab )), S_d), MATMUL(z_tld_k1_stab,m_x_k1_stab) )
    x_tld_k_of_zk1_stab   = MATMUL( MATMUL( TRANSPOSE(MATMUL(z_k1_stab , m_x_tld_k1_stab )), S_d), MATMUL(z_k1_stab,m_x_tld_k1_stab) )
      
    x_tld_k_of_xk1_dual     = MATMUL( MATMUL( m_x_k1_dual,     y_tld_k1_dual ) , MATMUL( z_tld_k1_dual , m_x_k1_dual     ) )
    x_tld_k_of_yk1_dual     = MATMUL( MATMUL( m_x_tld_k1_dual, y_k1_dual     ) , MATMUL( z_tld_k1_dual , m_x_tld_k1_dual ) )
    
    xp_tld_k_stab_gateaux = (x_tld_k_stab-x_tld_k_of_xk1_stab)/SQRT(SUM(dx_stab**2))
    zp_tld_k_stab_gateaux = (x_tld_k_stab-x_tld_k_of_zk1_stab)/SQRT(SUM(dz_stab**2))


    yp_tld_k_dual_gateaux = (x_tld_k_dual-x_tld_k_of_yk1_dual)/SQRT(SUM(dy_dual**2))







    FpX_stab  = MATMUL( TRANSPOSE(dz_k_dx_k1_stab) , MATMUL( s_d, z_tld_k_stab ) ) + MATMUL( TRANSPOSE(z_tld_k_stab) , MATMUL( s_d, dz_k_dx_k1_stab ) )

    FpZ_stab  = MATMUL( TRANSPOSE(dz_k_dz_k1_stab) , MATMUL( s_d, z_tld_k_stab ) ) + MATMUL( TRANSPOSE(z_tld_k_stab) , MATMUL( s_d, dz_k_dz_k1_stab ) )

    FpY_dual  = MATMUL( dy_k_dy_k1_dual , z_tld_k_dual    ) 

    ! dx~_k/dz_k-1 = ( < y~_k | z~_k > - < y~_k | z_k-1 >  m[ x~_k-1 ] )/dz_k-1
    zp_tld_k_dual_gateaux=(x_tld_k_dual-MATMUL( y_tld_k_dual,MATMUL(z_k1_dual,m_x_tld_k1_dual)))/SQRT(SUM(dz_dual**2))
    ! x~'_k_dz_k-1 =  < y~_k | dz^_k > m[ x~_k-1 ] 
    FpZ_dual  = MATMUL( y_tld_k_dual , MATMUL( dz_hat_dual , m_x_tld_k1_dual ) )
    ! dx~_k/dx_k-1 = ( < y~_k | z~_k >  - { m[x_k-1]<y~_k-1|z~_k> + <y~_k|z~_k-1>m[x_k-1] )/dx_k-1

    xp_tld_k_dual_gateaux=(x_tld_k_dual &
!                                      -(MATMUL(MATMUL(m_x_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
!                                       + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_k1_dual))) &
-      0.5d0*  (        MATMUL(MATMUL(m_x_tld_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
                + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_tld_k1_dual)) ) &
                            )/SQRT(SUM(dx_dual**2))

    WRITE(*,*)' x_tld_k_dual',SQRT(SUM(x_tld_k_dual**2))/  &
         SQRT(SUM(  (         MATMUL(MATMUL(m_x_tld_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
                            + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_tld_k1_dual)) )**2 )), & 
    SQRT(SUM(  (  MATMUL(MATMUL(m_x_tld_k1_dual,y_tld_k1_dual),z_tld_k_dual)  &
                + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_tld_k1_dual)) )**2 ))/SQRT(SUM( ( &
              MATMUL(MATMUL(m_x_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
            + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_k1_dual)))**2 )), &
            SQRT(SUM(x_tld_k_dual**2))/SQRT( SUM(( &
              MATMUL(MATMUL(m_x_k1_dual,y_tld_k1_dual),z_tld_k_dual))**2 &
            + SUM(MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_k1_dual)))**2 ))



    xp_tld_k_dual_gateaux=(  (MATMUL(MATMUL(m_x_tld_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
                            + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_tld_k1_dual))) &
                            -(MATMUL(MATMUL(m_x_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
                            + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_k1_dual))) &
                            )/SQRT(SUM(dx_dual**2))


    ! x~'_k_dx_k-1 = m[ x~_k-1 ] dx^_k-1 <y~_k-1 | dz~_k > + < y~_k | dz~_k-1 > m[ x~_k-1 ] dx^_k-1
    FpX_dual  = MATMUL( MATMUL( mp_dual , y_tld_k1_dual ), z_tld_k_dual    )  & 
              + MATMUL( y_tld_k_dual , MATMUL( z_tld_k1_dual , mp_dual) )

    FillN_stab=0d0
    FillN_dual=0d0
    do j=1,size(x_tld_k_stab,1)
       FillN_stab=FillN_stab+x_tld_k_stab(j,j)
       FillN_dual=FillN_dual+x_tld_k_dual(j,j)
    enddo
    FillN_stab=(FillN_stab-dble(n))/dble(n)
    FillN_dual=(FillN_dual-dble(n))/dble(n)

    WRITE(*,*)'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
    WRITE(*,*)' '
    WRITE(*,*)' (Tr<x_stab>-n)/n = ', FillN_stab
    WRITE(*,*)' '
    WRITE(*,*)'  dx_stab = ',sqrt(sum(dx_stab**2)),'  dz_stab = ',sqrt(sum(dz_stab**2))
    WRITE(*,*)' '
    WRITE(*,*)" ||f'_dz_k1_stab||   = ",SQRT(SUM(fpz_stab**2)), SQRT(SUM( zp_tld_k_stab_gateaux**2))
    WRITE(*,*)" ||f'_dx_k1_stab||   = ",SQRT(SUM(fpx_stab**2)), SQRT(SUM( xp_tld_k_stab_gateaux**2))
    WRITE(*,*)' '
    WRITE(*,*)' (Tr<x_dual>-n)/n = ',FillN_dual
    WRITE(*,*)' '
    WRITE(*,*)'  dx_dual = ',sqrt(sum(dx_dual**2)),'  dz_dual = ',sqrt(sum(dz_dual**2)),'  dy_dual = ',sqrt(sum(dy_dual**2))
    WRITE(*,*)' '
    WRITE(*,*)" ||f'_dy_k1_dual||   = ",SQRT(SUM(fpy_dual**2)), SQRT(SUM( yp_tld_k_dual_gateaux**2)) 
    WRITE(*,*)" ||f'_dz_k1_dual||   = ",SQRT(SUM(fpz_dual**2)), SQRT(SUM( zp_tld_k_dual_gateaux**2)) 
    WRITE(*,*)" ||f'_dx_k1_dual||   = ",SQRT(SUM(fpx_dual**2)),SQRT(SUM( xp_tld_k_dual_gateaux**2)) 
    WRITE(*,*)' '
    WRITE(*,*)'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '

    WRITE(99,44)kount, tau, FillN_stab, sqrt(sum(dx_stab**2)),sqrt(sum(dz_stab**2)),SQRT(SUM(fpz_stab**2)),SQRT(SUM(fpx_stab**2)), &
                            FillN_dual, sqrt(sum(dx_dual**2)),sqrt(sum(dz_dual**2)),sqrt(sum(dy_dual**2)),SQRT(SUM(fpy_dual**2)), &
                                        SQRT(SUM(fpz_dual**2)),SQRT(SUM(fpx_dual**2))

44  format(I3,"   ", 20(e12.6,",   "))

    deallocate(dz_stab)
    deallocate(dz_hat_stab)
    deallocate(dx_stab)
    deallocate(dx_hat_stab)
    deallocate(dy_dual)     
    
    deallocate(dy_hat_dual) 
    deallocate(dz_dual)     
    deallocate(dz_hat_dual) 
    deallocate(dx_dual)     
    deallocate(dx_hat_dual) 
    deallocate(mp_stab) 
    deallocate(mp_dual) 
    deallocate(dz_k_dx_k1_stab)
    deallocate(dz_k_dz_k1_dual)
    deallocate(dy_k_dy_k1_dual)
    deallocate(dz_k_dx_k1_dual)
    deallocate(dy_k_dx_k1_dual)
    deallocate(x_tld_k_of_xk1_stab)
    deallocate(x_tld_k_of_zk1_stab)      
    deallocate(x_tld_k_of_xk1_dual)
    deallocate(x_tld_k_of_yk1_dual)
    deallocate(x_tld_k_of_zk1_dual)   
    deallocate(xp_tld_k_stab_gateaux) 
    deallocate(zp_tld_k_stab_gateaux) 
    deallocate(xp_tld_k_dual_gateaux) 
    deallocate(yp_tld_k_dual_gateaux) 
    deallocate(zp_tld_k_dual_gateaux) 
    deallocate(FpX_stab)  
    deallocate(FpZ_stab)  
    deallocate(FpX_dual)  
    deallocate(FpZ_dual)  
    deallocate(FpY_dual)  

  END SUBROUTINE SpAMMsand_Error_Analysis
#endif

  FUNCTION SpAMMsand_Basis_Compare(a,b) RESULT(compare)
    REAL(SpAMM_kind)                   :: compare, MAX_DOT
    REAL(spamm_kind), dimension(:,:)   :: a,b
    compare=SQRT( SUM( ( MATMUL(a,b) - MATMUL(b,a) )**2 ) ) & 
           /SQRT( SUM(a**2)) /SQRT( SUM(b**2) )
  END FUNCTION SpAMMsand_Basis_Compare

!#ifdef DENSE_DIAGNOSTICS

  FUNCTION spammsand_shift_tree_2d( x, low_prev, high_prev, low_new, high_new, stab ) RESULT(d)

!#else
!
!  FUNCTION spammsand_shift_tree_2d( x, low_prev, high_prev, low_new, high_new ) RESULT(d)
!
!#endif

    !!!!!    shft=low_new+(x-low_prev)*(high_new-low_new)/(high_prev-low_prev)

    TYPE(spamm_tree_2d_symm) ,   POINTER     :: d
    TYPE(spamm_tree_2d_symm) ,   POINTER     :: x
    REAL(SpAMM_KIND), OPTIONAL,INTENT(IN)    :: low_prev, high_prev, low_new, high_new   
    REAL(SpAMM_KIND)                         :: SHFT,SCAL 

!#ifdef DENSE_DIAGNOSTICS
    logical :: stab
!#endif

    integer :: i

    SHFT=low_new-low_prev*(high_new-low_new)/(high_prev-low_prev)
    SCAL=(high_new-low_new)/(high_prev-low_prev)

#ifdef DENSE_DIAGNOSTICS
    IF(stab)then
       scal_shift_stab=scal
       shft_shift_stab=shft       
       m_x_k1_stab=m_x_k1_stab*scal
       m_x_k1_stab=m_x_k1_stab*shft*i_d
    ELSE
       scal_shift_dual=scal
       shft_shift_dual=shft       
       m_x_k1_dual=m_x_k1_dual*scal
       m_x_k1_dual=m_x_k1_dual+shft*i_d
    ENDIF
#endif

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

  END FUNCTION spammsand_shift_tree_2d

#ifdef DENSE_DIAGNOSTICS

  FUNCTION spammsand_scaled_invsqrt_mapping( x, sc, stab ) result(d) 

#else

  FUNCTION spammsand_scaled_invsqrt_mapping( x, sc ) result(d) 

#endif

    TYPE(spamm_tree_2d_symm), POINTER  :: d
    TYPE(spamm_tree_2d_symm), POINTER  :: x
    REAL(SpAMM_KIND),      INTENT(IN)  :: sc 
    REAL(SpAMM_KIND)                   :: SHFT,SCAL 

#ifdef DENSE_DIAGNOSTICS
    logical :: stab
#endif

    SHFT=SpAMM_half*SQRT(sc)*SpAMM_three
    SCAL=SpAMM_half*(-sc)*SQRT(sc)

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

#ifdef DENSE_DIAGNOSTICS
    IF(stab)then
       scal_mapp_stab=scal
       shft_mapp_stab=shft       
       m_x_k1_stab=m_x_k1_stab*scal
       m_x_k1_stab=m_x_k1_stab+shft*i_d
       m_x_tld_k1_stab=m_x_tld_k1_stab*scal
       m_x_tld_k1_stab=m_x_tld_k1_stab+shft*i_d

!       WRITE(*,*)' stab scal shift = ',scal,shft
!       WRITE(*,*)'stab d = ',d%frill%norm2, SUM(m_x_k1_stab**2)
    ELSE
       scal_mapp_dual=scal
       shft_mapp_dual=shft       
       m_x_k1_dual=m_x_k1_dual*scal       
       m_x_k1_dual=m_x_k1_dual+shft*i_d
       m_x_tld_k1_dual=m_x_tld_k1_dual*scal       
       m_x_tld_k1_dual=m_x_tld_k1_dual+shft*i_d

!       WRITE(*,*)' dual scal shift = ',scal,shft
!       WRITE(*,*)'dual d = ',d%frill%norm2, SUM(m_x_k1_dual**2)
    ENDIF
#endif

  END FUNCTION spammsand_scaled_invsqrt_mapping

  FUNCTION spammsand_scaling_invsqrt(xo) RESULT(sc)

    REAL(SpAMM_KIND) :: xo, sc
    sc=MIN( Approx3, SpAMM_three/( SpAMM_one + SQRT(xo) + xo) )    

  END FUNCTION spammsand_scaling_invsqrt
 
end module SpAMMsand_inverse_squareroot

! Nested SpAMM solvers (SpAMM sandwitches) for matrix functions
program SpAMM_sandwich_inverse_squareroot

  USE spammpack 
  USE sandpack 
  USE SpAMMsand_inverse_squareroot
  use test_utilities
  
  implicit none

  real(spamm_kind), dimension(:,:), ALLOCATABLE :: s_dense

  TYPE(spammsand_tree_2d_slices), pointer        :: z, z_head
!, y, y_head

  type(SpAMM_tree_2d_symm),       pointer        :: s => null()
  type(SpAMM_tree_2d_symm),       pointer        :: x => null()
  type(SpAMM_tree_2d_symm),       pointer        :: t => null()
  type(SpAMM_tree_2d_symm),       pointer        :: s_orgnl => null()
  type(SpAMM_tree_2d_symm),       pointer        :: z_total => null()

!  real(spamm_kind), dimension(:, :), allocatable :: S_dense

  character(len = 1000)                          :: matrix_filename
  real(SpAMM_KIND)                               :: x_hi, x_new, logtau_strt, logtau_stop, logtau_dlta, & 
                                                    tau_dlta, tau_xtra, error, tmp1,tmp2, final_tau, s_work, zs_work
  logical :: first, second

  integer, parameter                             :: slices=6

  real(SpAMM_KIND), dimension(1:slices)          :: tau
 
  integer :: i,j,k, kount

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)
  S_dense=SpAMM_half*(S_dense+TRANSPOSE(S_dense))


  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm( S_DENSE, in_O = s )

#ifdef DENSE_DIAGNOSTICS
  N = s%frill%ndimn(1)
  LWORK = 1+6*N+2*N**2
  LIWORK = 3+5*N    
  allocate(eval(N))
  allocate(work(LWORK))
  allocate(iwork(LIWORK))
#endif
  !=============================================================
  ! the max eigenvalue

  x_hi = SpAMMSand_rqi_extremal(s,1d-4,high_O=.TRUE.)
  WRITE(*,*)' hi extremal = ',x_hi

  ! normalize the max ev of s to 1.  
  s       => SpAMM_scalar_times_tree_2d_symm( SpAMM_one/x_hi , s )
  s_orgnl => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, in_O = s_orgnl, threshold_O = SpAMM_normclean )

  logtau_strt=-2                                      ! starting accuracy
  logtau_stop=-10                                      ! stoping  "
  logtau_dlta=(logtau_stop-logtau_strt)/dble(slices-1) ! span (breadth) of SpAMM thresholds 
  tau_dlta=10d0**logtau_dlta
  final_tau=10d0**logtau_stop                          ! penultimate spamm threshold
    
  allocate(z) 
  z_head => z ! top of the sandwich
  do i=1,slices

     z%tau = 10d0**( logtau_strt + logtau_dlta * float(i-1) )
     tau(i)=z%tau
     z%mtx => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )

     if(i==slices)then        
        z%nxt => null()
     else
        allocate(z%nxt) 
        z => z%nxt
     endif

  enddo

   write(*,33)tau
33 format(' building |Z> = ',4('|',e6.1,'>'),'...')

  ! work matrices ...
  x => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
  t => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )

  kount=1
!  first=.FALSE.
  first=.TRUE.
  second=.TRUE.

  z=>z_head

  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>
     
     ! WATER block=32, tau=0.07
     ! BAD TUBE block=32, tau=0.01

     z%tau = 0.01d0
     call spammsand_scaled_newton_shulz_inverse_squareroot( s, x, z%mtx, z%tau, first, second, kount )

     STOP

     first=.TRUE.
     second=.FALSE.
     if(.not.associated(z%nxt))exit    

     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z_total, z%mtx, z%nxt%tau , NT_O=.TRUE., in_O = t )
     z_total => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = z_total, threshold_O = z%nxt%tau )      

     ! Nested sandwich, to recompute the error every time. 
     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z_total, s_orgnl, z%nxt%tau*1d-2, NT_O=.FALSE., in_O = t )
     s => SpAMM_tree_2d_symm_times_tree_2d_symm(       t, z_total, z%nxt%tau     , NT_O=.TRUE. , in_O = s )

      s_work=s%frill%flops/dble(s%frill%ndimn(1))**3
     zs_work=t%frill%flops/dble(t%frill%ndimn(1))**3

     write(*,*)' t work = ',zs_work,', s_work = ',s_work

     x_new =  SpAMMSand_rqi_extremal( x, z%nxt%tau*1d-2 , high_O=.TRUE. )
     WRITE(*,*)' hi extremal = ',x_new 

     x     => SpAMM_scalar_times_tree_2d_symm( SpAMM_one / x_new , x )
     x_hi  = x_hi * x_new

     s => SpAMM_tree_2d_symm_copy_tree_2d_symm( x, in_O = s , threshold_O = SpAMM_normclean )

     z => z%nxt
     
  enddo

  ! de-normalize, back to max_ev (1->max_ev)...  
  z=>z_head
  z%mtx => SpAMM_scalar_times_tree_2d_symm(SpAMM_one/SQRT(x_hi), z%mtx)

  ! at this point, we are basically done.  the rest is IO/verification of trace, error, timing & etc
  ! timers down ...

  ! x <= |z_spammsand> = |z_1> . |z_2> ... |z_slices>.  Has to be applied as left (T) and right (N) (its not symmetric)
  x => SpAMM_tree_2d_symm_copy_tree_2d_symm( z%mtx, in_O = x )
  z => z%nxt
  do while(associated(z)) ! build the inverse factors |z> = |z_1>.|z_2> ... |z_slices> (right handed)
     t => SpAMM_tree_2d_symm_times_tree_2d_symm( x, z%mtx, 1d-20, nt_O=.TRUE., in_O = t )     
     x => SpAMM_tree_2d_symm_copy_tree_2d_symm( t , in_O = x, threshold_O = SpAMM_zero )
     z => z%nxt
  enddo

  ! get the original back ...
  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 
  
  ! I = <z^T_spammsand|s|z_spammsand>
  t => SpAMM_tree_2d_symm_times_tree_2d_symm( s, x, 1d-20, nt_O=.TRUE. , in_O = t )   
  s => SpAMM_tree_2d_symm_times_tree_2d_symm( x, t, 1d-20, nt_O=.FALSE., in_O = s )   

  ! error = [ <Zt_1|Zt_2|...<Zt_slice| S |Z_slice>...|Z_2>|Z_1> -1.*I ] /N**2

  s => SpAMM_scalar_plus_tree_2d_symm( -SpAMM_one, s) 

  error=SQRT(s%frill%norm2)/dble(s%frill%ndimn(1))**2

  write(*,*)' error = ',error


end program SpAMM_sandwich_inverse_squareroot


























!!$
!!$  FUNCTION SpAMMsand_zee_dot_x_error( z, x, tau) RESULT(d)
!!$
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN) :: z, x
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER             :: d
!!$    REAL(spamm_kind) :: tau 
!!$              
!!$    real(spamm_kind), dimension(1:z%frill%ndimn(1),1:z%frill%ndimn(2)) :: x_d, z_d, d_d
!!$    
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_d )
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_d )
!!$
!!$    d => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau , nt_O=.TRUE., in_O = in_o )
!!$
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( d, d_d )
!!$
!!$    d_d = d_d - MATMUL( z_d , x_d )
!!$
!!$!    d_d = 0.5d0*(d_d+TRANSPOSE(d_d))
!!$
!!$    d => SpAMM_convert_dense_to_tree_2d_symm( d_d, in_o = d )
!!$
!!$  END FUNCTION SpAMMsand_zee_dot_x_error

!!$
!!$  FUNCTION SpAMMsand_compare_naiv(z,s,x) RESULT(compare)
!!$
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN) :: z, s, x
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER             :: d
!!$
!!$    real(spamm_kind)                               :: compare           
!!$
!!$    real(spamm_kind), dimension(1:s%frill%ndimn(1),1:s%frill%ndimn(2)) :: s_d, z_d, x_naiv, x_d
!!$    
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d )
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_d )
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_d )
!!$
!!$    x_naiv=MATMUL(z_d,MATMUL(s_d,z_d))
!!$    x_d=x_d-x_naiv
!!$
!!$    compare=SpAMMsand_Basis_Compare_d(s_d,x_d)
!!$
!!$  END FUNCTION SpAMMsand_compare_naiv
!!$
!!$  FUNCTION SpAMMsand_compare_stab(z,s,x) RESULT(compare)
!!$
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN) :: z, s,x
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER             :: d
!!$
!!$    real(spamm_kind)                               :: compare           
!!$
!!$    real(spamm_kind), dimension(1:s%frill%ndimn(1),1:s%frill%ndimn(2)) :: s_d, z_d, x_naiv, x_d
!!$    
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d )
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_d )
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_d )
!!$
!!$    x_naiv=MATMUL(TRANSPOSE(z_d),MATMUL(s_d,z_d))
!!$    x_d=x_d-x_naiv
!!$
!!$    compare=SpAMMsand_Basis_Compare_d(s_d,x_d)
!!$
!!$  END FUNCTION SpAMMsand_compare_stab
!!$



!!$
!!$
!!$    MAX_DOT=-1D10
!!$    DO I=1,N
!!$       DO J=1,N
!!$          MAX_DOT=MAX(MAX_DOT,DOT_PRODUCT(A_D(:,I),B_D(:,J)))
!!$       ENDDO
!!$    ENDDO
!!$
!!$    WRITE(*,*)'a',MAX_DOT,ACOS(MAX_DOT)
!!$
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( a, a_d )
!!$    a_d=TRANSPOSE(a_d)
!!$    call dsyevd("V", "U", N, a_d, N, eval, work, LWORK, iwork, LIWORK, info)
!!$
!!$    MAX_DOT=-1D10
!!$    DO I=1,N
!!$       DO J=1,N
!!$          MAX_DOT=MAX(MAX_DOT,DOT_PRODUCT(A_D(:,I),B_D(:,J)))
!!$       ENDDO
!!$    ENDDO
!!$
!!$    WRITE(*,*)'b', MAX_DOT,ACOS(MAX_DOT)
!!$
!!$
!!$
!!$!    svd=MATMUL(TRANSPOSE(a_d),b_d)
!!$!    call dsyevd("V", "U", N, svd, N, eval, work, LWORK, iwork, LIWORK, info)
!!$!    WRITE(*,*)' ev = ',eval(1),acos(eval(1))
!!$
!!$!    svd=MATMUL(a_d,b_d)
!!$!    call dsyevd("V", "U", N, svd, N, eval, work, LWORK, iwork, LIWORK, info)
!!$!    WRITE(*,*)' ev = ',eval(1:2)
!!$
!!$
!!$!A = orth(A);
!!$!B = orth(B);
!!$!s = svd(A'*B);
!!$!a = acos(min(s,1));
!!$
!!$    MAX_DOT=-1D10
!!$    DO I=1,N
!!$       DO J=1,N
!!$          MAX_DOT=MAX(MAX_DOT,DOT_PRODUCT(A_D(:,I),B_D(:,J)))
!!$       ENDDO
!!$    ENDDO
!!$
!!$    WRITE(*,*)MAX_DOT,ACOS(MAX_DOT)
!!$

!!$    compare =ACOS(MAX_DOT)



!!$
!!$  FUNCTION SpAMMsand_dense_stab(z,s, in_o) RESULT(d)
!!$
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN) :: z, s
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER, OPTIONAL   :: in_o
!!$    TYPE(SpAMM_tree_2d_symm) , POINTER             :: d
!!$               
!!$    real(spamm_kind), dimension(1:s%frill%ndimn(1),1:s%frill%ndimn(2)) :: s_d, z_d, x_stab
!!$    
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d )
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_d )
!!$
!!$    x_stab=MATMUL(TRANSPOSE(z_d),MATMUL(s_d,z_d))
!!$    d=> SpAMM_convert_dense_to_tree_2d_symm( x_stab, in_o = in_o )
!!$
!!$  END FUNCTION SpAMMsand_dense_stab

