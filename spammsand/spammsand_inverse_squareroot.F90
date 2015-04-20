#define DENSE_DIAGNOSTICS
! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan  -Wextra -std=f2008 "

module SpAMMsand_inverse_squareroot
  USE spammpack 

  implicit none

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1

#ifdef DENSE_DIAGNOSTICS
  real(spamm_kind), dimension(:,:), ALLOCATABLE :: &
     i_d,s_d,z_k,z_k1,z_tld_k,z_tld_k1,  y_k,y_k1,y_tld_k,y_tld_k1,x_k,x_k1,x_tld_k_naiv, &
     x_tld_k_stab,x_tld_k1,m_x_k1,m_x_tld_k1,x_tld_k_yz 

  REAL(spamm_kind) :: scal_shift,shft_shift, scal_mapp, shft_mapp
  
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

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, t, tau, DoDuals, second, kount)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER :: x, z, y, t
!    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: xx, z, t
    REAL(SpAMM_KIND) :: Tau
!    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau
!    LOGICAL, INTENT(IN)                              :: first
    LOGICAL                                           :: DoDuals,second
    INTEGER, INTENT(INOUT)                            :: kount
    INTEGER                                           :: i,  j, k
    REAL(SpAMM_KIND)                                  :: sc, TrX, tau_xtra
    REAL(SpAMM_KIND)                                  :: x_work,z_work,y_work,zs_work,sz_work
    REAL(SpAMM_KIND)                                  :: x_fill,z_fill,y_fill,zs_fill,sz_fill

    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev
    REAL(SpAMM_KIND)                                  :: tau_zdotz,sz_norm

#ifdef DENSE_DIAGNOSTICS
#else
    IF(DoDuals)tHeN
#endif
       ! y_0 => s
       y => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       y => SpAMM_tree_2d_symm_copy_tree_2d_symm( S , in_o = y, threshold_O = SpAMM_normclean ) 
#ifdef DENSE_DIAGNOSTICS
#else
    endIF
#endif
    ! z_0 => I
    z=>SpAMM_set_identity_2d_symm ( s%frill%ndimn, in_o = z )
    !  x_0 => s
    x => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_O = x, threshold_O = SpAMM_normclean ) 

#ifdef DENSE_DIAGNOSTICS
    n=x%frill%ndimn(1)
    ALLOCATE(i_d     (1:N,1:N))
    ALLOCATE(s_d     (1:N,1:N))
    ALLOCATE(z_k     (1:N,1:N))
    ALLOCATE(z_k1    (1:N,1:N))
    ALLOCATE(z_tld_k (1:N,1:N))
    ALLOCATE(z_tld_k1(1:N,1:N))
    ALLOCATE(y_k     (1:N,1:N))
    ALLOCATE(y_k1    (1:N,1:N))
    ALLOCATE(y_tld_k (1:N,1:N))
    ALLOCATE(y_tld_k1(1:N,1:N))
    ALLOCATE(m_x_k1  (1:N,1:N))
    ALLOCATE(m_x_tld_k1  (1:N,1:N))
    ALLOCATE(x_k     (1:N,1:N))
    ALLOCATE(x_k1    (1:N,1:N))
    ALLOCATE(x_tld_k1(1:N,1:N))
    ALLOCATE(x_tld_k_naiv (1:N,1:N))
    ALLOCATE(x_tld_k_stab (1:N,1:N))
    ALLOCATE(x_tld_k_yz (1:N,1:N))
    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d  )
    CALL SpAMM_convert_tree_2d_symm_to_dense( z, i_d  )
    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_k1 )
    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_tld_k1 )
    CALL SpAMM_convert_tree_2d_symm_to_dense( s, y_k1 )
    CALL SpAMM_convert_tree_2d_symm_to_dense( s, y_tld_k1 )
    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_k1 )
    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k1 )
    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k_stab )
    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k_naiv )
    CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k_yz )
    DoDuals=.TRUE.
#endif
    kount=1
    sz_norm=1d0
    FillN=1d10
    DO i = 0, 26

       ! check the trace for convergence:
       FillN_prev=FillN
       TrX=SpAMM_trace_tree_2d_symm_recur(x)
       FillN = abs( dble(s%frill%ndimn(1)) - TrX )/dble(s%frill%ndimn(1))       

       WRITE(*,33)tau, kount, TrX, FillN, 1d2*x%frill%non0s/dble(x%frill%ndimn(1))**2, y_work*1d2 , z_work*1d2 , sz_work*1d2 , x_work*1d2 
33     format('  ... Tr< ',e6.1,', i=',i2,' > = ', F18.10,' dN=',e10.3,' %of N^2 = ',f10.5,'%,  %of N^3 = ',f10.5,'%, ',f10.5,'%, ',f10.5,'%, ',f10.5,'%')

       IF(i>2 .and. FillN<0.1d0 .AND. FillN>FillN_prev )then

!          WRITE(*,*)' fill n = ',filln,' filln_prev ',filln_prev
!          write(*,*)' elevation' 

          RETURN  ! Elevation

       endif
       IF( FillN <  Tau )then

!          WRITE(*,*)' fill n = ',filln,' tau**2 = ',tau
!          write(*,*)' anhiliaiton '

          RETURN  ! Anihilation

       end IF

#ifdef DENSE_DIAGNOSTICS
       m_x_k1=x_k1
       scal_shift=1d0
       shft_shift=0d0 
       scal_mapp =1d0
       shft_mapp =0d0
#endif

       IF(FillN>0.4d0)THEN
          delta=10.d-2 ! maybe this should be a variable too, passed in?
          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_zero)
       ELSEIF(FillN>0.1d0)THEN
          delta=1.0d-2 ! maybe this should be a variable too, passed in?
          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_half)
       ELSE
          sc=1d0
       ENDIF


       x => spammsand_scaled_invsqrt_mapping( x, sc )

       tau_zdotz=1d-3

!tau*SQRT(SQRT(z%frill%norm2))
!       WRITE(*,*)' tauz = ',tau_zdotz


       ! |Z_n+1> =  <Z_n| X_n>  
       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau_zdotz , nt_O=.TRUE., in_O = t )
       ! update (copy+threshold)
       z => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = z, threshold_O = tau_zdotz) 
       ! stats
       z_work=t%frill%flops/dble(t%frill%ndimn(1))**3
       z_fill=z%frill%non0s

#ifdef DENSE_DIAGNOSTICS
       CALL SpAMM_convert_tree_2d_symm_to_dense( x, m_x_tld_k1 )
       CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_tld_k )
       z_k = MATMUL( z_k1, m_x_k1 )
#else
       IF(DoDuals)tHeN
#endif
          ! <Y_n+1|
          t => SpAMM_tree_2d_symm_times_tree_2d_symm( x, y, tau , nt_O=.TRUE., in_O = t )
          y_work=t%frill%flops/dble(t%frill%ndimn(1))**3

          ! update (+threshold)
          y => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = y, threshold_O = tau ) 

          ! <X_n> = <Y_n|Z_n>
          x => SpAMM_tree_2d_symm_times_tree_2d_symm( y, z, tau   , nt_O=.TRUE. , in_O = x )
          x_work=x%frill%flops/dble(x%frill%ndimn(1))**3

#ifdef DENSE_DIAGNOSTICS
          y_k = MATMUL( m_x_k1,  y_k1 )
          CALL SpAMM_convert_tree_2d_symm_to_dense( y, y_tld_k )
          CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k_yz )
#else
       ELSE
#endif
          !   |t> =      S|Z_k> 
          t => SpAMM_tree_2d_symm_times_tree_2d_symm( s, z,  tau*1d-4  , NT_O=.TRUE. , in_O = t )
          sz_work=t%frill%flops/dble(t%frill%ndimn(1))**3
          sz_norm=SQRT(t%frill%norm2)
          
#ifdef DENSE_DIAGNOSTICS
          x => SpAMM_tree_2d_symm_times_tree_2d_symm( z, t,  tau  , NT_O=.TRUE. , in_O = x )
          CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k_naiv )
#endif

          ! <X_k> = <Z_k|S|Z_k>
          x => SpAMM_tree_2d_symm_times_tree_2d_symm( z, t,  tau      , NT_O=.FALSE. , in_O = x )
          x_work=x%frill%flops/dble(x%frill%ndimn(1))**3
          
#ifdef DENSE_DIAGNOSTICS
          CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_tld_k_stab )          
          x_k = MATMUL( TRANSPOSE( z_k ) , MATMUL( s_d, z_k ) )

          if(i>0)CALL SpAMMsand_Error_Analysis(kount, tau, FillN)
          y_k1=y_k
          z_k1=z_k
          x_k1=x_k
          y_tld_k1=Y_tld_k 
          Z_tld_k1=Z_tld_k 
          x_tld_k1=x_tld_k_stab 
#else
       endif
#endif
       kount=kount+1
    END DO

    CaLl SpAMM_destruct_tree_2d_symm_recuR (Y)
   
  END SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot

#ifdef DENSE_DIAGNOSTICS
  SUBROUTINE SpAMMsand_Error_Analysis(kount, tau, FillN)

    integer :: i,kount

    real(spamm_kind), dimension(:,:), ALLOCATABLE :: &
         MP, MP_Gateaux, DX,dx_hat, zp_k,yp_k, fp_naiv, fp_stab, fp_yz, &
         x_tld_k_of_xk1_naiv, x_tld_k_of_xk1_stab, x_tld_k_of_xk1_yz, &
         xp_tld_k_naiv_gateaux,xp_tld_k_stab_gateaux,xp_tld_k_yz_gateaux, dz, dz_hat

    real(spamm_kind) :: x_tld_k_naiv_compare, x_tld_k_stab_compare, x_tld_k_natv_compare, tau, filln

    ALLOCATE(dX  (1:N,1:N))
    ALLOCATE(dX_hat(1:N,1:N))

    ALLOCATE(dZ  (1:N,1:N))
    ALLOCATE(dZ_hat(1:N,1:N))

    ALLOCATE(Mp  (1:N,1:N))
    ALLOCATE(Mp_Gateaux  (1:N,1:N))
    ALLOCATE(zp_k(1:N,1:N))
    ALLOCATE(yp_k(1:N,1:N))
    ALLOCATE(fp_naiv      (1:N,1:N))
    ALLOCATE(fp_stab      (1:N,1:N))
    ALLOCATE(fp_yz      (1:N,1:N))
    ALLOCATE( x_tld_k_of_xk1_naiv (1:N,1:N))
    ALLOCATE( x_tld_k_of_xk1_stab (1:N,1:N))
    ALLOCATE( x_tld_k_of_xk1_yz (1:N,1:N))
    ALLOCATE(xp_tld_k_naiv_gateaux (1:N,1:N))
    ALLOCATE(xp_tld_k_stab_gateaux (1:N,1:N))
    ALLOCATE(xp_tld_k_yz_gateaux (1:N,1:N))

    dx        = x_tld_k1 - x_k1
    dx_hat    = dx/SQRT(SUM(dx**2))

    dz        = z_tld_k1 - z_k1
    dz_hat    = dz/SQRT(SUM(dz**2))

    mp        = scal_mapp * scal_shift * dx_hat
    mp_Gateaux=(m_x_tld_k1-m_x_k1)/SQRT(SUM(dx**2))

    zp_k = MATMUL( z_tld_k1, mp ) 
    yp_k = MATMUL( mp, y_tld_k1 ) 

    x_tld_k_of_xk1_yz     = MATMUL( MATMUL( m_x_k1, y_tld_k1 ) , &
                            MATMUL( z_tld_k1 , m_x_k1 ) )
    x_tld_k_of_xk1_naiv   = MATMUL( MATMUL(           MATMUL(z_tld_k1 , m_x_k1 ) , S_d), &
                            MATMUL(z_tld_k1,m_x_k1) )
    x_tld_k_of_xk1_stab   = MATMUL( MATMUL( TRANSPOSE(MATMUL(z_tld_k1 , m_x_k1 )), S_d), &
                            MATMUL(z_tld_k1,m_x_k1) )
    
    xp_tld_k_yz_gateaux   = (x_tld_k_yz  -x_tld_k_of_xk1_yz  )/SQRT(SUM(dx**2))
    xp_tld_k_naiv_gateaux = (x_tld_k_naiv-x_tld_k_of_xk1_naiv)/SQRT(SUM(dx**2))
    xp_tld_k_stab_gateaux = (x_tld_k_stab-x_tld_k_of_xk1_stab)/SQRT(SUM(dx**2))

    fp_yz    = MATMUL( yp_k , z_tld_k ) + MATMUL( y_tld_k , zp_k )
    fp_naiv  = MATMUL(     zp_k, MATMUL( s_d,  z_tld_k ) ) + MATMUL(  z_tld_k, MATMUL( s_d,     zp_k ) )
    fp_stab  = MATMUL(TRANSPOSE(zp_k),MATMUL(s_d,z_tld_k))+MATMUL(TRANSPOSE(z_tld_k),MATMUL(s_d,zp_k))

!    WRITE(*,*)" ||f'_yz     ||   = ",SQRT(SUM(fp_yz**2)),   SQRT(SUM( xp_tld_k_yz_gateaux**2)) 
!    WRITE(*,*)" ||f'_naiv   ||   = ",SQRT(SUM(fp_naiv**2)), SQRT(SUM( xp_tld_k_naiv_gateaux**2)) 
!    WRITE(*,*)" ||f'_stab   ||   = ",SQRT(SUM(fp_stab**2)), SQRT(SUM( xp_tld_k_stab_gateaux**2))

!    WRITE(*,*)" ||f'_naiv   ||   = ",SQRT(SUM(fp_naiv**2)), SQRT(SUM( xp_tld_k_naiv_gateaux**2)) 
!    WRITE(*,*)" ||f'_stab   ||   = ",SQRT(SUM(fp_stab**2)), SQRT(SUM( xp_tld_k_stab_gateaux**2))

!    IF(kount==2)WRITE(*,43)
!    WRITE(*,44)kount, tau, FillN, SQRT(SUM(dX**2))                   , &
!         SQRT(SUM(fp_yz  **2)), SQRT(SUM( xp_tld_k_yz_gateaux  **2)) , &
!         SQRT(SUM(fp_naiv**2)), SQRT(SUM( xp_tld_k_naiv_gateaux**2)) , &
!         SQRT(SUM(fp_stab**2)), SQRT(SUM( xp_tld_k_stab_gateaux**2))
 

    fp_stab  = MATMUL(TRANSPOSE(dz),MATMUL(s_d,z_tld_k))+MATMUL(TRANSPOSE(z_tld_k),MATMUL(s_d,dz))

    WRITE(*,*)" ||f'_stab   ||   = ",SQRT(SUM(dz**2)),SQRT(SUM(fp_stab**2))

















   
43     FORMAT(" Itr          Tau,             %N,            /\X,        f'_stab,       /\f_stab,        f'_dual,       /\f_dual,        f'_naiv,       /\f_naiv ")
44     format(I3,"   ", 20(e12.6,",   "))

    deALLOCATE(dX        )
    deALLOCATE(dX_hat    )

    deALLOCATE(dZ        )
    deALLOCATE(dZ_hat    )

    deALLOCATE(Mp        )
    deALLOCATE(Mp_Gateaux)
    deALLOCATE(zp_k      )
    deALLOCATE(yp_k      )
    deALLOCATE(fp_naiv   )
    deALLOCATE(fp_stab   )
    deALLOCATE(fp_yz     )
    deALLOCATE( x_tld_k_of_xk1_naiv  )
    deALLOCATE( x_tld_k_of_xk1_stab  )
    deALLOCATE( x_tld_k_of_xk1_yz  )
    deALLOCATE(xp_tld_k_naiv_gateaux )
    deALLOCATE(xp_tld_k_stab_gateaux )
    deALLOCATE(xp_tld_k_yz_gateaux )

  END SUBROUTINE SpAMMsand_Error_Analysis
#endif

  FUNCTION SpAMMsand_Basis_Compare(a,b) RESULT(compare)
    REAL(SpAMM_kind)                   :: compare, MAX_DOT
    REAL(spamm_kind), dimension(:,:)   :: a,b
    compare=SQRT( SUM( ( MATMUL(a,b) - MATMUL(b,a) )**2 ) ) & 
           /SQRT( SUM(a**2)) /SQRT( SUM(b**2) )
  END FUNCTION SpAMMsand_Basis_Compare

  FUNCTION spammsand_shift_tree_2d( x, low_prev, high_prev, low_new, high_new ) RESULT(d)
    !!!!!    shft=low_new+(x-low_prev)*(high_new-low_new)/(high_prev-low_prev)

    TYPE(spamm_tree_2d_symm) ,   POINTER     :: d
    TYPE(spamm_tree_2d_symm) ,   POINTER     :: x
    REAL(SpAMM_KIND), OPTIONAL,INTENT(IN)    :: low_prev, high_prev, low_new, high_new   
    REAL(SpAMM_KIND)                         :: SHFT,SCAL 

    integer :: i

    SHFT=low_new-low_prev*(high_new-low_new)/(high_prev-low_prev)
    SCAL=(high_new-low_new)/(high_prev-low_prev)

#ifdef DENSE_DIAGNOSTICS
    scal_shift=scal
    shft_shift=shft       
    m_x_k1=m_x_k1*scal
    m_x_k1=m_x_k1+shft*i_d
#endif

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

  END FUNCTION spammsand_shift_tree_2d

  FUNCTION spammsand_scaled_invsqrt_mapping( x, sc ) result(d) 
    TYPE(spamm_tree_2d_symm), POINTER  :: d
    TYPE(spamm_tree_2d_symm), POINTER  :: x
    REAL(SpAMM_KIND),      INTENT(IN)  :: sc 
    REAL(SpAMM_KIND)                   :: SHFT,SCAL 

    SHFT=SpAMM_half*SQRT(sc)*SpAMM_three
    SCAL=SpAMM_half*(-sc)*SQRT(sc)

#ifdef DENSE_DIAGNOSTICS
    scal_mapp=scal
    shft_mapp=shft
    m_x_k1=m_x_k1*scal
    m_x_k1=m_x_k1+shft*i_d
#endif

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

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

  integer, parameter                             :: slices=4

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

  logtau_strt=-9                                       ! starting accuracy
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
  first=.FALSE.
  second=.TRUE.

  z=>z_head
  ! z_total => I
  z_total=>SpAMM_set_identity_2d_symm ( s%frill%ndimn, in_o = z_total )

  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>
     
     call spammsand_scaled_newton_shulz_inverse_squareroot( s, x, z%mtx, t, z%tau, first, second, kount )

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

