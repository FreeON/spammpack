!#define DENSE_DIAGNOSTICS
! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan  -Wextra -std=f2008 "

module SpAMMsand_inverse_squareroot

  USE spammpack 

  implicit none

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1
  !
#ifdef DENSE_DIAGNOSTICS


  character(len=132) :: file_dual, file_stab, corename

  real(spamm_kind), dimension(:,:), ALLOCATABLE ::   i_d, s_d, y_k, y_k1, z_k, z_k1, x_k, x_k1, m_x_k1,          &   
       z_tld_k_stab, z_tld_k1_stab, x_tld_k_stab, x_tld_k1_stab,  m_x_tld_k1_stab,  y_tld_k_stab, y_tld_k1_stab, &
       z_tld_k_dual, z_tld_k1_dual, x_tld_k_dual, x_tld_k1_dual,  m_x_tld_k1_dual,  y_tld_k_dual, y_tld_k1_dual

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

  FUNCTION DblToChar(D)
    CHARACTER(LEN=*),PARAMETER  :: INTERNAL_DBL_FMT='(D10.4)'
    REAL*8,INTENT(IN)         :: D
    CHARACTER(LEN=22) :: DblToChar
    WRITE(UNIT=DblToChar,FMT=INTERNAL_DBL_FMT)D
    DblToChar=ADJUSTL(DblToChar)
  END FUNCTION DblToChar

  FUNCTION IntToChar(I)
    CHARACTER(LEN=*), PARAMETER :: INTERNAL_INT_FMT='(I22)'
    INTEGER  :: I
    CHARACTER(LEN=22) :: IntToChar
    WRITE(UNIT=IntToChar,FMT=INTERNAL_INT_FMT)I
    IntToChar=ADJUSTL(IntToChar)
  END FUNCTION IntToChar

  FUNCTION CharToInt(C)
    CHARACTER(LEN=*), PARAMETER :: INTERNAL_INT_FMT='(I22)'
    CHARACTER(LEN=*),INTENT(IN) :: C
    INTEGER                     :: CharToInt
    READ(UNIT=C,FMT=INTERNAL_INT_FMT,ERR=666) CharToInt
    RETURN
666 WRITE(*,*)"[CharToInt] Fatal error, can not convert "//TRIM(C)
  END FUNCTION CharToInt

  SUBROUTINE LowCase(String)
    CHARACTER(LEN=*), PARAMETER :: Upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(LEN=*), PARAMETER :: Lower='abcdefghijklmnopqrstuvwxyz'
    INTEGER                     :: I,J
    CHARACTER(LEN=*)            :: String
    DO I=1,LEN(String)
      J=INDEX(Upper,String(I:I))
      IF(J.NE.0)String(I:I)=Lower(J:J)
    ENDDO
  END SUBROUTINE LowCase

  FUNCTION CharToDbl(C)
    CHARACTER(LEN=*),INTENT(IN) :: C
    REAL*8                      :: CharToDbl
    INTEGER                     :: TempInt
    CHARACTER(LEN=*),PARAMETER  :: INTERNAL_DBL_FMT='(D22.16)'
!!$    CALL LowCase(C)
!!$    IF(SCAN(C,'.') == 0 .AND. SCAN(C,'e') == 0 .AND. SCAN(C,'d') == 0)THEN
!!$      READ(UNIT=C,FMT=INTERNAL_INT_FMT,ERR=667)TempInt
!!$      CharToDbl=DBLE(TempInt)
!!$    ELSEIF(SCAN(C,".") == 0 .AND. (SCAN(C,'e') /= 0 .OR. SCAN(C,'d') /= 0)) THEN
!!$      STOP "Check floating point number input format for "//TRIM(C)
!!$    ELSE
      READ(UNIT=C,FMT=INTERNAL_DBL_FMT,ERR=667)CharToDbl

!!$    ENDIF
    RETURN
667 STOP ' Oooops in chartodbl '
  END FUNCTION CharToDbl

  FUNCTION Sigmoid(Scale, Inflect, x)
    REAL(SpAMM_KIND) :: Scale, Inflect, Sigmoid, x
    Sigmoid=SpAMM_One/(SpAMM_One+EXP(- Scale * (x-Inflect) ) )

  END FUNCTION Sigmoid

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, Tau_0, Tau_S, delta_0, &
                                                              DoDuals, RightTight, DoScale, First, kount)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN) :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER             :: x,z ! OUT
    TYPE(SpAMM_tree_2d_symm) , POINTER             :: x_stab, z_stab, y_stab, x_dual, z_dual, y_dual, y_tmp, z_tmp
    REAL(SpAMM_KIND)                               :: Tau_0, Tau_S, delta_0
    LOGICAL                                        :: DoDuals, DoScale, First, RightTight
    INTEGER                                        :: kount
    INTEGER                                        :: i,  j, k
    REAL(SpAMM_KIND)                               :: scale, delta, TrX, FillN, FillN_prev
    REAL(SpAMM_KIND)                               :: y_stab_work,z_stab_work,x_stab_work, y_dual_work,z_dual_work,x_dual_work
    REAL(SpAMM_KIND)                               :: y_stab_fill,z_stab_fill,x_stab_fill, y_dual_fill,z_dual_fill,x_dual_fill
    CHARACTER                                      :: RT

#ifdef DENSE_DIAGNOSTICS
    INTEGER                                        :: stat
    !
    !
    file_dual=TRIM(corename)//'_t0'//TRIM(DblToChar(tau_0))//'_ts'//TRIM(DblToChar(tau_s)) &
            //'_d'//TRIM(DblToChar(delta_0))//'_b'//TRIM(IntToChar(SpAMM_BLOCK_SIZE))//'_dual.dat'

    file_stab=TRIM(corename)//'_t0'//TRIM(DblToChar(tau_0))//'_ts'//TRIM(DblToChar(tau_s)) &
            //'_d'//TRIM(DblToChar(delta_0))//'_b'//TRIM(IntToChar(SpAMM_BLOCK_SIZE))//'_'//TRIM(RT)//'_stab.dat'
    !
    open(unit=98, iostat=stat, file=file_stab,status='old')
    if(stat.eq.0) close(98, status='delete')
    open(unit=98, iostat=stat, file=file_stab,status='new')
    close(98)

    open(unit=99, iostat=stat, file=file_dual,status='old')
    if(stat.eq.0) close(99, status='delete')
    open(unit=99, iostat=stat, file=file_dual,status='new')
    close(99)   

#else
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

#ifdef DENSE_DIAGNOSTICS
#else
    ELSE
#endif
       ! y_0 => s
       y_stab => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       y_stab => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_o = y_stab, threshold_O = SpAMM_normclean ) 
       !  z_0 => I
       z_stab => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       z_stab => SpAMM_set_identity_2d_symm( s%frill%ndimn, in_o = z_stab )
       !  x_0 => s
       x_stab => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       x_stab => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_O = x_stab, threshold_O = SpAMM_normclean ) 

#ifdef DENSE_DIAGNOSTICS
#else
    endIF
#endif

    y_tmp  => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
    z_tmp  => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )

#ifdef DENSE_DIAGNOSTICS
    n=s%frill%ndimn(1)
    ALLOCATE(i_d(1:N,1:N));  ALLOCATE(s_d(1:N,1:N))
    allocate(y_k(1:n,1:n));  allocate(y_k1(1:n,1:n));  
    allocate(z_k(1:n,1:n));  allocate(z_k1(1:n,1:n));  
    allocate(x_k(1:n,1:n));  allocate(x_k1(1:n,1:n));  
    allocate(m_x_k1(1:n,1:n))
    allocate(z_tld_k_stab (1:n,1:n)); allocate(z_tld_k1_stab  (1:n,1:n));  allocate(x_tld_k_stab(1:n,1:n)); 
    allocate(x_tld_k1_stab(1:n,1:n)); allocate(m_x_tld_k1_stab(1:n,1:n)); allocate(y_tld_k_stab(1:n,1:n));   allocate(y_tld_k1_stab(1:n,1:n))         
    allocate(z_tld_k_dual (1:n,1:n)); allocate(z_tld_k1_dual  (1:n,1:n));  allocate(x_tld_k_dual(1:n,1:n));  
    allocate(x_tld_k1_dual(1:n,1:n)); allocate(m_x_tld_k1_dual (1:n,1:n)); allocate(y_tld_k_dual(1:n,1:n));  allocate(y_tld_k1_dual(1:n,1:n)) 

    ! I=diag(1)
    i_d=SpAMM_zero
    DO k=1,n
       i_d(k,k)=SpAMM_one
    ENDDO
    !
    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d  )

    y_k=s_d; y_k1=s_d; 
    x_k=s_d; x_k1=s_d
    z_k=i_d; z_k1=i_d; 
    y_tld_k1_stab=s_d; z_tld_k1_stab=i_d; x_tld_k1_stab=s_d
    y_tld_k1_dual=s_d; z_tld_k1_dual=i_d; x_tld_k1_dual=s_d

#endif

    kount=1
    FillN=1d10
    DO i = 0, 28

#ifdef DENSE_DIAGNOSTICS
       IF(I>0)THEN
          TrX=SpAMM_trace_tree_2d_symm_recur(x_dual)
          FillN = MIN(FillN, abs( dble(s%frill%ndimn(1)) - TrX )/dble(s%frill%ndimn(1)))       
          WRITE(*,33)tau_0, delta, scale, kount, TrX, FillN, y_dual_work*1d2 , z_dual_work*1d2 , x_dual_work*1d2           
          TrX=SpAMM_trace_tree_2d_symm_recur(x_stab)
          FillN = MIN(FillN, abs( dble(s%frill%ndimn(1)) - TrX )/dble(s%frill%ndimn(1)))       
          WRITE(*,34)tau_0, delta, scale, RightTight, kount, TrX, FillN, y_stab_work*1d2 , z_stab_work*1d2 , x_stab_work*1d2 
       ENDIF
       FillN=0d0
       do j=1,n
          FillN=FillN+x_k(j,j)
       enddo
       FillN=(dble(n)-FillN)/dble(n)
#else
       ! check the trace for convergence:
       FillN_prev=FillN
       IF(DoDuals)THEN
          TrX=SpAMM_trace_tree_2d_symm_recur(x_dual)
       ELSE
          TrX=SpAMM_trace_tree_2d_symm_recur(x_stab)
       ENDIF

       FillN = MIN(FillN, abs( dble(s%frill%ndimn(1)) - TrX )/dble(s%frill%ndimn(1)))       

       IF(I>0)THEN
          IF(DoDuals)THEN
             WRITE(*,33)tau_0, delta, scale, kount, TrX, FillN, y_dual_work*1d2 , z_dual_work*1d2 , x_dual_work*1d2 
          elsE
             WRITE(*,34)tau_0, delta, scale, RightTight, kount, TrX, FillN, y_stab_work*1d2 , z_stab_work*1d2 , x_stab_work*1d2 
          ENDIF
       ENDIF
#endif

33     format('  dual Tr<t0=',e8.3,',d=',e8.3,',s=',e8.3,',    n=',i2,' > = ', F18.10,' dN=',e10.3, & 
              ', y_wrk: ',f10.5,'%, z_wrk:',f10.5,'%, x_wrk:',f10.5,'%')
34     format('  stab Tr<t0=',e8.3,',d=',e8.3,',s=',e8.3,',r=',L1,',n=',i2,' > = ', F18.10,' dN=',e10.3, &
              ', y_wrk: ',f10.5,'%, z_wrk:',f10.5,'%, x_wrk:',f10.5,'%')

       ! convergence
       IF(i>2 .and. FillN<0.1d0 .AND. FillN>FillN_prev )then
          !          WRITE(*,*)' fill n = ',filln,' filln_prev ',filln_prev
          !          write(*,*)' elevation' 
          !          RETURN  ! Elevation
       endif

       IF( FillN <  Tau_0 )then
          !          WRITE(*,*)' fill n = ',filln,' tau_0**2 = ',tau
          !          write(*,*)' anhiliaiton '
          !          RETURN  ! Anihilation
       end IF

       ! stabilization (only on first call)
       IF(First)THEN
          delta=delta_0*Sigmoid(75d0, 3.d-1, FillN)
       ELSE
          delta=SpAMM_Zero
       ENDIF
! .35
!   dual Tr<t0=.100E-01,d=.758E-16,s=.100E+01,    n=15 > =    3100.2924960903 dN= 0.195E-04, y_wrk:    5.98005%, z_wrk:   4.49350%, x_wrk:   2.99249%
! .3
!  dual Tr<t0=.100E-01,d=.113E-13,s=.100E+01,    n=15 > =    3100.2811656597 dN= 0.237E-04, y_wrk:    5.94334%, z_wrk:   4.44845%, x_wrk:   2.99481%

       
       ! scaling 
       IF(DoScale)THEN
          scale=SpAMM_One + ( Approx3 - SpAMM_One )*Sigmoid(50d0, 3.5d-1, FillN)
       ELSE       
          scale=SpAMM_One
       ENDIF

#ifdef DENSE_DIAGNOSTICS       
       scal_shift=1d0
       shft_shift=0d0 
       scal_mapp =1d0
       shft_mapp =0d0
#else
       IF(DoDuals)THEN
#endif
          ! m[x_n-1,c]   
          x_dual => spammsand_shift_tree_2d( x_dual, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          x_dual => spammsand_scaled_invsqrt_mapping( x_dual, scale)
          ! |z_n> =  |z_n-1> m[x_n-1]                    
          z_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_dual, x_dual, tau_0, nt_O=.TRUE., & 
                   in_O = z_tmp , stream_file_O='z_dual'//inttoCHAR(I) )
          z_dual=> SpAMM_tree_2d_symm_copy_tree_2d_symm( z_tmp, in_O = z_dual, threshold_O = tau_0 ) 
          z_dual_work=z_tmp%frill%flops/dble(z_tmp%frill%ndimn(1))**3; z_dual_fill=z%frill%non0s
#ifdef DENSE_DIAGNOSTICS       
#else
       ELSE 
#endif
          ! m[x_n-1,c]   
          x_stab => spammsand_shift_tree_2d( x_stab, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          x_stab => spammsand_scaled_invsqrt_mapping( x_stab, scale)
          ! |z_n> =  |z_n-1> m[x_n-1]  
          z_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( z_stab, x_stab, tau_0, nt_O=.TRUE., in_O = z_tmp , &
                   stream_file_O='z_stab'//inttoCHAR(I) )
          z_stab=> SpAMM_tree_2d_symm_copy_tree_2d_symm(  z_tmp, in_O = z_stab, threshold_O = tau_0 ) 
          z_stab_work=z_tmp%frill%flops/dble(z_tmp%frill%ndimn(1))**3
          z_stab_fill=z%frill%non0s
#ifdef DENSE_DIAGNOSTICS       
#else
       ENDIF
#endif

#ifdef DENSE_DIAGNOSTICS
       CALL SpAMM_convert_tree_2d_symm_to_dense( z_stab, z_tld_k_stab )
       CALL SpAMM_convert_tree_2d_symm_to_dense( z_dual, z_tld_k_dual )
#else
       IF(DoDuals)tHeN
#endif
          ! <y_n| = m[x_n-1]<y_n-1|
          y_tmp  => SpAMM_tree_2d_symm_times_tree_2d_symm( x_dual, y_dual, tau_S , nt_O=.TRUE. , &
                    in_O = y_tmp , stream_file_O='y_dual_'//inttoCHAR(I) )
          y_dual => SpAMM_tree_2d_symm_copy_tree_2d_symm( y_tmp, in_O = y_dual, threshold_O = tau_S )
          ! x_n = <y_n|z_n>
          x_dual => SpAMM_tree_2d_symm_times_tree_2d_symm( y_dual, z_dual, tau_0 , nt_O=.TRUE. , &
                    in_O = x_dual )
          ! stats ...   ! double chck that stats for ytmp==ydual
          y_dual_work=y_tmp%frill%flops/dble(y_tmp%frill%ndimn(1))**3   ; y_dual_fill=y_tmp%frill%non0s
          x_dual_work=x_dual%frill%flops/dble(x_dual%frill%ndimn(1))**3 ; x_dual_fill=x_dual%frill%non0s
#ifdef DENSE_DIAGNOSTICS
          CALL SpAMM_convert_tree_2d_symm_to_dense( y_dual, y_tld_k_dual )
          CALL SpAMM_convert_tree_2d_symm_to_dense( x_dual, x_tld_k_dual )
#else
       ELSE
#endif
          IF(RightTight)THEN
             RT='R' ! R flag
             ! | w_n > = < s | z_n >
             y_stab => SpAMM_tree_2d_symm_times_tree_2d_symm( s     , z_stab, tau_S , NT_O=.TRUE.  ,   &
                       in_O = y_stab )
             ! | x_n > = < zt_n | w_n >
             x_stab => SpAMM_tree_2d_symm_times_tree_2d_symm( z_stab, y_stab, tau_0 , NT_O=.FALSE. ,   &
                       in_O = x_stab , stream_file_O='x_stab_'//inttoCHAR(I) )
             ! stats ...
             y_stab_work=y_stab%frill%flops/dble(y_stab%frill%ndimn(1))**3;  y_stab_fill=y_stab%frill%non0s
             x_stab_work=x_stab%frill%flops/dble(x_stab%frill%ndimn(1))**3;  x_stab_fill=x_stab%frill%non0s
          ELSE
             RT='L' ! L flag
             ! | y_n > = < zt_n | s >
             y_stab => SpAMM_tree_2d_symm_times_tree_2d_symm( z_stab, s     , tau_S , NT_O=.FALSE. ,   &
                       in_O = y_stab )
             ! | x_n > = < y_n | z_n >
             x_stab => SpAMM_tree_2d_symm_times_tree_2d_symm( y_stab, z_stab, tau_0 , NT_O=.TRUE.  ,   &
                       in_O = x_stab , stream_file_O='x_stab_'//inttoCHAR(I) )
             ! stats ...
             y_stab_work=y_stab%frill%flops/dble(y_stab%frill%ndimn(1))**3; y_stab_fill=y_stab%frill%non0s
             x_stab_work=x_stab%frill%flops/dble(x_stab%frill%ndimn(1))**3; x_stab_fill=x_stab%frill%non0s
          ENDIF
#ifdef DENSE_DIAGNOSTICS
          CALL SpAMM_convert_tree_2d_symm_to_dense( y_stab, y_tld_k_stab )
          CALL SpAMM_convert_tree_2d_symm_to_dense( x_stab, x_tld_k_stab )          
          calL SpAMMsand_Error_Analysis(kount, tau_0, tau_S, scale, delta)
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

  SUBROUTINE SpAMMsand_Error_Analysis(kount, tau, tau_xtra, sc, delta)

     real(spamm_kind) :: Tau, Tau_xtra, sc, delta, FillN_stab, FillN_dual, &
         dy_stab_sclr,dz_stab_sclr, &
         dx_stab_sclr,dy_dual_sclr, &
         dz_dual_sclr,dx_dual_sclr

    integer :: i,kount,j,k

    real(spamm_kind), dimension(:,:), ALLOCATABLE ::                        &
         dy_stab,dy_hat_stab,dz_stab,dz_hat_stab,dx_stab,dx_hat_stab,       &
         dy_dual,dy_hat_dual,dz_dual,dz_hat_dual,dx_dual,dx_hat_dual,       &
         mp_dy_stab,mp_dz_stab,mp_dy_dual,mp_dz_dual,mp_dx_stab,mp_dx_dual, &
         dz_dy_stab,dy_dy_stab,dx_dy_stab,dy_dz_stab,dz_dz_stab,dx_dz_stab, &
         dz_dy_dual,dy_dy_dual,dx_dy_dual,dy_dz_dual,dz_dz_dual,dx_dz_dual, &
         dx_dx_stab, dx_dx_dual, x_k_of_z, x_k_of_y

    ALLOCATE(dy_stab(1:N,1:N) );   ALLOCATE(dy_hat_stab(1:N,1:N) )
    ALLOCATE(dz_stab(1:N,1:N) );   ALLOCATE(dz_hat_stab(1:N,1:N) )
    ALLOCATE(dx_stab(1:N,1:N) );   ALLOCATE(dx_hat_stab(1:N,1:N) )
    ALLOCATE(dy_dual(1:N,1:N) );   ALLOCATE(dy_hat_dual(1:N,1:N) )
    ALLOCATE(dz_dual(1:N,1:N) );   ALLOCATE(dz_hat_dual(1:N,1:N) )
    ALLOCATE(dx_dual(1:N,1:N) );   ALLOCATE(dx_hat_dual(1:N,1:N) )
    ALLOCATE(mp_dy_stab(1:N,1:N) );ALLOCATE(mp_dz_stab(1:N,1:N) ); 
    ALLOCATE(mp_dy_dual(1:N,1:N) );ALLOCATE(mp_dz_dual(1:N,1:N) ); 
    ALLOCATE(dz_dy_stab(1:N,1:N) );ALLOCATE(dy_dy_stab(1:N,1:N) )
    ALLOCATE(dx_dy_stab(1:N,1:N) );ALLOCATE(dy_dz_stab(1:N,1:N) )
    ALLOCATE(dz_dz_stab(1:N,1:N) );ALLOCATE(dx_dz_stab(1:N,1:N) )
    ALLOCATE(dz_dy_dual(1:N,1:N) );ALLOCATE(dy_dy_dual(1:N,1:N) )
    ALLOCATE(dx_dy_dual(1:N,1:N) );ALLOCATE(dy_dz_dual(1:N,1:N) )
    ALLOCATE(dz_dz_dual(1:N,1:N) );ALLOCATE(dx_dz_dual(1:N,1:N) )
    ALLOCATE(dx_dx_stab(1:N,1:N) );ALLOCATE(dx_dx_dual(1:N,1:N) )
    ALLOCATE(mp_dx_stab(1:N,1:N) );ALLOCATE(mp_dx_dual(1:N,1:N) )
    ALLOCATE(x_k_of_y(1:N,1:N) );  ALLOCATE(x_k_of_z(1:N,1:N) );

    !-----------------
    dy_stab     = y_tld_k1_stab - y_k1;   dy_hat_stab = dy_stab/SQRT(SUM(dy_stab**2))
    dz_stab     = z_tld_k1_stab - z_k1;   dz_hat_stab = dz_stab/SQRT(SUM(dz_stab**2))
    dx_stab     = x_tld_k1_stab - x_k1;   dx_hat_stab = dx_stab/SQRT(SUM(dx_stab**2))    
    dy_stab_sclr=SQRT(SUM(dy_stab**2))
    dz_stab_sclr=SQRT(SUM(dz_stab**2))
    dx_stab_sclr=SQRT(SUM(dx_stab**2))
    ! ------
    dy_dual     = y_tld_k1_dual - y_k1;   dy_hat_dual = dy_dual/SQRT(SUM(dy_dual**2))
    dz_dual     = z_tld_k1_dual - z_k1;   dz_hat_dual = dz_dual/SQRT(SUM(dz_dual**2))
    dx_dual     = x_tld_k1_dual - x_k1;   dx_hat_dual = dx_dual/SQRT(SUM(dx_dual**2))
    dy_dual_sclr=SQRT(SUM(dy_dual**2))
    dz_dual_sclr=SQRT(SUM(dz_dual**2))
    dx_dual_sclr=SQRT(SUM(dx_dual**2))
    !---------------------------------------------
    ! unperturbed ...
    m_x_k1=x_k1
    ! apply the shift and scale by delta (stabilization) ...
    m_x_k1=m_x_k1*scal_shift
    m_x_k1=m_x_k1+shft_shift*i_d
    ! apply the NS map ...
    m_x_k1=m_x_k1*scal_mapp
    m_x_k1=m_x_k1+shft_mapp*i_d
    ! here is the reference dual, y & k ...
    y_k = MATMUL(m_x_k1,y_k1)
    z_k = MATMUL(z_k1,m_x_k1)
    ! taking the "dual" reference in double precision ... 
    x_k = MATMUL(y_k,    z_k)
    !   x_k = MATMUL(MATMUL(TRANSPOSE(z_k),s_d),z_k)
    !---------------------------------------------
    mp_dy_stab  = scal_mapp * scal_shift * MATMUL(dy_hat_stab,z_k1)
    mp_dz_stab  = scal_mapp * scal_shift * MATMUL(y_k1,dz_hat_stab)
    mp_dx_stab  = scal_mapp * scal_shift * dx_hat_stab
    ! ------
    mp_dy_dual  = scal_mapp * scal_shift * MATMUL(dy_hat_dual,z_k1)
    mp_dz_dual  = scal_mapp * scal_shift * MATMUL(y_k1,dz_hat_dual)
    mp_dx_dual  = scal_mapp * scal_shift * dx_hat_stab

    ! ------
    mp_dx_stab  = scal_mapp * scal_shift * dx_hat_stab

    mp_dz_stab  = scal_mapp * scal_shift * ( MATMUL(TRANSPOSE(dz_hat_stab),MATMUL(S_d,Z_k)) &
                                            +MATMUL(TRANSPOSE(z_k),MATMUL(S_d,dZ_hat_stab)) )

    dz_dz_stab   = MATMUL( dz_hat_stab , m_x_k1 ) + MATMUL( z_k1 , mp_dz_stab )

    dx_dz_stab   = MATMUL(TRANSPOSE(dz_dz_stab),MATMUL(s_d,z_k))+&
                   MATMUL(TRANSPOSE(z_k),MATMUL(s_d,dz_dz_stab))

    ! ------
    dy_dy_dual   = MATMUL(m_x_k1,dy_hat_dual)+MATMUL(mp_dy_dual,y_k1)
    dz_dy_dual   = MATMUL(z_k1,mp_dy_dual)
    dx_dy_dual   = MATMUL(dy_dy_dual,z_k) + MATMUL(y_k,dz_dy_dual)

    ! ------
    dy_dz_dual   = MATMUL(mp_dz_dual,y_k1)
    dz_dz_dual   = MATMUL(dz_hat_dual,m_x_k1)+MATMUL(z_k1,mp_dz_dual)
    dx_dz_dual   = MATMUL(dy_dz_dual,z_k) + MATMUL(y_k,dz_dz_dual)

    ! ------
    dx_dx_stab   = MATMUL(mp_dx_stab,MATMUL(y_k1,z_k))+MATMUL(y_k,MATMUL(z_k1,mp_dx_stab))
    dx_dx_dual   = MATMUL(mp_dx_dual,MATMUL(y_k1,z_k))+MATMUL(y_k,MATMUL(z_k1,mp_dx_dual))

    FillN_stab=0d0
    FillN_dual=0d0
    do j=1,size(x_tld_k_stab,1)
       FillN_stab=FillN_stab+x_tld_k_stab(j,j)
       FillN_dual=FillN_dual+x_tld_k_dual(j,j)
    enddo
    FillN_stab=abs(dble(n)-FillN_stab)/dble(n)
    FillN_dual=abs(dble(n)-FillN_dual)/dble(n)

    WRITE(*,*)'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
    WRITE(*,*)'   stab: '
    WRITE(*,23)SQRT(SUM(dx_dz_stab**2)), SQRT(SUM(dz_stab**2)),SQRT(SUM(dx_dz_stab**2))*SQRT(SUM(dz_stab**2))
    WRITE(*,21)SQRT(SUM(dx_dx_stab**2)), SQRT(SUM(dx_stab**2)),SQRT(SUM(dx_dx_stab**2))*SQRT(SUM(dx_stab**2))
    WRITE(*,24)FillN_stab
    WRITE(*,*)'   dual: '
    WRITE(*,22)SQRT(SUM(dx_dy_dual**2)), SQRT(SUM(dy_dual**2)),SQRT(SUM(dx_dy_dual**2))*SQRT(SUM(dy_dual**2))
    WRITE(*,23)SQRT(SUM(dx_dz_dual**2)), SQRT(SUM(dz_dual**2)),SQRT(SUM(dx_dz_dual**2))*SQRT(SUM(dz_dual**2))
    WRITE(*,21)SQRT(SUM(dx_dx_dual**2)), SQRT(SUM(dx_dual**2)),SQRT(SUM(dx_dx_dual**2))*SQRT(SUM(dx_dual**2))
    WRITE(*,24)FillN_dual

22  FORMAT("    f'_dy_k1 = ",e12.6,", dy_k1 = ",e12.6,", f'_dy*dy_k1 = ",e12.6)
23  FORMAT("    f'_dz_k1 = ",e12.6,", dz_k1 = ",e12.6,", f'_dz*dz_k1 = ",e12.6)
21  FORMAT("    f'_dx_k1 = ",e12.6,", dx_k1 = ",e12.6,", f'_dx*dx_k1 = ",e12.6)
24  FORMAT("    [n-trx]/n = ",e12.6)
    !
    open(unit=98, file=file_stab, position='append')
    WRITE(98,44)kount, tau, tau_xtra, sc, delta, FillN_stab, &
         0.d0, 0d0, &
         SQRT(SUM(dx_dz_stab**2)), SQRT(SUM(dz_stab**2)),  &
         SQRT(SUM(dx_dx_stab**2)), SQRT(SUM(dx_stab**2))
    close(98)

    open(unit=99, file=file_dual, position='append')
    WRITE(99,44)kount, tau, tau_xtra, sc, delta, FillN_dual, &
    SQRT(SUM(dx_dy_dual**2)), SQRT(SUM(dy_dual**2)), SQRT(SUM(dx_dz_dual**2)), SQRT(SUM(dz_dual**2)),  &
         SQRT(SUM(dx_dx_dual**2)), SQRT(SUM(dx_dual**2))

    close(99)

44  format(I3,"   ", 20(e12.6,",   "))

    ! free some memory ...
    DEALLOCATE(x_k_of_y);   DEALLOCATE(x_k_of_z); 
    DEALLOCATE(dy_stab);    DEALLOCATE(dy_hat_stab); 
    DEALLOCATE(dz_stab);    DEALLOCATE(dz_hat_stab)
    DEALLOCATE(dx_stab);    DEALLOCATE(dx_hat_stab)
    DEALLOCATE(dy_dual);    DEALLOCATE(dy_hat_dual)
    DEALLOCATE(dz_dual);    DEALLOCATE(dz_hat_dual)
    DEALLOCATE(dx_dual);    DEALLOCATE(dx_hat_dual)
    DEALLOCATE(mp_dy_stab); DEALLOCATE(mp_dz_stab);  DEALLOCATE(mp_dx_stab)
    DEALLOCATE(mp_dy_dual); DEALLOCATE(mp_dz_dual);  DEALLOCATE(mp_dx_dual)
    DEALLOCATE(dz_dy_stab); DEALLOCATE(dy_dy_stab)
    DEALLOCATE(dx_dy_stab); DEALLOCATE(dy_dz_stab)
    DEALLOCATE(dz_dz_stab); DEALLOCATE(dx_dz_stab)
    DEALLOCATE(dz_dy_dual); DEALLOCATE(dy_dy_dual)
    DEALLOCATE(dx_dy_dual); DEALLOCATE(dy_dz_dual)
    DEALLOCATE(dz_dz_dual); DEALLOCATE(dx_dz_dual)
    DEALLOCATE(dx_dx_stab); DEALLOCATE(dx_dx_dual)

    ! update the privous iterates ...
    y_k1=y_k
    z_k1=z_k
    x_k1=x_k
    y_tld_k1_dual=y_tld_k_dual
    z_tld_k1_dual=z_tld_k_dual
    x_tld_k1_dual=x_tld_k_dual
    y_tld_k1_stab=y_tld_k_stab
    z_tld_k1_stab=z_tld_k_stab
    x_tld_k1_stab=x_tld_k_stab
    !

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

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

#ifdef DENSE_DIAGNOSTICS
    scal_mapp=scal
    shft_mapp=shft       
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

  logical :: DoDuals, DoScale, First, RightTight
  real(SpAMM_KIND)                               :: tau_0, tau_S, delta

  integer, parameter                             :: slices=6

  real(SpAMM_KIND), dimension(1:slices)          :: tau

  integer :: i,j,k, kount
  character(len=10) :: c_tau_0, c_tau_S,  c_scale, c_delta, c_dual, c_righttight


  !  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)
  call get_command_argument(2, c_tau_0)
  call get_command_argument(3, c_tau_S)
  call get_command_argument(4, c_delta)
  call get_command_argument(5, c_dual)
  call get_command_argument(6, c_scale)
  call get_command_argument(7, c_righttight)

  tau_0=CharToDbl(c_tau_0)
  tau_S=CharToDbl(c_tau_S)
  delta=CharToDbl(c_delta)
  if(ADJUSTL(c_righttight)=='R')then
     RightTight=.TRUE.
  else
     RightTight=.FALSE.
  endif

  call read_MM(matrix_filename, S_dense)
  S_dense=SpAMM_half*(S_dense+TRANSPOSE(S_dense))

  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm( S_DENSE, in_O = s )



#ifdef DENSE_DIAGNOSTICS
  call get_command_argument(8, corename)

  N = s%frill%ndimn(1)
  LWORK = 1+6*N+2*N**2
  LIWORK = 3+5*N    
  allocate(eval(N))
  allocate(work(LWORK))
  allocate(iwork(LIWORK))
#endif

  s => SpAMM_scalar_plus_tree_2d_symm(  1d-2, s)
  !=============================================================
  ! the max eigenvalue

  x_hi = SpAMMSand_rqi_extremal(s,1d-10,high_O=.TRUE.)



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

     z%tau_0 = tau_0 !10d0**( logtau_strt + logtau_dlta * float(i-1) )
     z%tau_S = tau_S !z%Tau_0*1d-2
     tau(i)=z%tau_0
     z%mtx => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )

     if(i==slices)then       
        z%nxt => null()
     else
        allocate(z%nxt) 
        z => z%nxt
     endif
  enddo

!  write(*,33)tau!!$
33 format(' building |Z> = ',4('|',e6.1,'>'),'...')

  ! work matrices ...
  x => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
  t => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )

  kount=1

  First=.TRUE.

  IF(c_dual=='D')THEN
     DoScale=.TRUE.
  ELSE
     DoDuals=.FALSE.
  ENDIF

  IF(c_scale=='S')THEN
     DoScale=.TRUE.
  ELSE
     DoScale=.FALSE.
  ENDIF

  z=>z_head

  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>

     ! For now, just hardening the preconditioner

     call spammsand_scaled_newton_shulz_inverse_squareroot( s, x, z%mtx, z%tau_0, z%tau_S, delta,  &
                                                            DoDuals, RightTight, DoScale, First, kount)

     STOP
!!$
!!$     first=.TRUE.
!!$     second=.FALSE.
!!$     if(.not.associated(z%nxt))exit    
!!$
!!$     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z_total, z%mtx, z%nxt%tau , NT_O=.TRUE., in_O = t )
!!$     z_total => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = z_total, threshold_O = z%nxt%tau )      
!!$
!!$     ! Nested sandwich, to recompute the error every time. 
!!$     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z_total, s_orgnl, z%nxt%tau*1d-2, NT_O=.FALSE., in_O = t )
!!$     s => SpAMM_tree_2d_symm_times_tree_2d_symm(       t, z_total, z%nxt%tau     , NT_O=.TRUE. , in_O = s )
!!$
!!$     s_work=s%frill%flops/dble(s%frill%ndimn(1))**3
!!$     zs_work=t%frill%flops/dble(t%frill%ndimn(1))**3
!!$
!!$     write(*,*)' t work = ',zs_work,', s_work = ',s_work
!!$
!!$     x_new =  SpAMMSand_rqi_extremal( x, z%nxt%tau*1d-2 , high_O=.TRUE. )
!!$     WRITE(*,*)' hi extremal = ',x_new 
!!$
!!$     x     => SpAMM_scalar_times_tree_2d_symm( SpAMM_one / x_new , x )
!!$     x_hi  = x_hi * x_new
!!$
!!$     s => SpAMM_tree_2d_symm_copy_tree_2d_symm( x, in_O = s , threshold_O = SpAMM_normclean )
!!$
!!$     z => z%nxt

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



!!$
!!$  SUBROUTINE SpAMMsand_Error_Analysis2(kount, tau)
!!$
!!$    real(spamm_kind) :: Tau, FillN_stab, FillN_dual
!!$
!!$    integer :: i,kount,j,k
!!$
!!$    real(spamm_kind), dimension(:,:), ALLOCATABLE :: &
!!$    dz_stab, dz_hat_stab, dx_stab, dx_hat_stab,       &
!!$    dy_dual, dy_hat_dual,dz_dual, dz_hat_dual,       &
!!$    dx_dual, dx_hat_dual,mp_stab, mp_dual,           &
!!$    dz_k_dx_k1_stab,dz_k_dz_k1_stab,dz_k_dz_k1_dual, &
!!$    dy_k_dy_k1_dual,dz_k_dx_k1_dual,dy_k_dx_k1_dual, &
!!$    x_tld_k_of_xk1_stab,x_tld_k_of_zk1_stab,         &
!!$    x_tld_k_of_xk1_dual,x_tld_k_of_yk1_dual,x_tld_k_of_zk1_dual,       &
!!$    xp_tld_k_stab_gateaux,zp_tld_k_stab_gateaux,xp_tld_k_dual_gateaux, & 
!!$    yp_tld_k_dual_gateaux,zp_tld_k_dual_gateaux,                       &
!!$    FpX_stab,FpZ_stab,FpX_dual,FpZ_dual,FpY_dual  
!!$
!!$    allocate(dz_stab                (1:N,1:N) )
!!$    allocate(dz_hat_stab            (1:N,1:N) )
!!$    allocate(dx_stab                (1:N,1:N) )
!!$    allocate(dx_hat_stab            (1:N,1:N) )
!!$
!!$    allocate(dy_dual                (1:N,1:N) )
!!$    allocate(dy_hat_dual            (1:N,1:N) )
!!$    allocate(dz_dual                (1:N,1:N) )
!!$    allocate(dz_hat_dual            (1:N,1:N) )
!!$    allocate(dx_dual                (1:N,1:N) )
!!$    allocate(dx_hat_dual            (1:N,1:N) )
!!$    allocate(mp_stab                (1:N,1:N) )
!!$    allocate(mp_dual                (1:N,1:N) )
!!$
!!$    allocate(dz_k_dx_k1_stab        (1:N,1:N) )
!!$    allocate(dz_k_dz_k1_stab        (1:N,1:N) )
!!$
!!$    allocate(dz_k_dz_k1_dual        (1:N,1:N) )
!!$    allocate(dy_k_dy_k1_dual        (1:N,1:N) )
!!$
!!$    allocate(dz_k_dx_k1_dual        (1:N,1:N) )
!!$    allocate(dy_k_dx_k1_dual        (1:N,1:N) )
!!$
!!$    allocate(x_tld_k_of_xk1_stab    (1:N,1:N) )
!!$    allocate(x_tld_k_of_zk1_stab    (1:N,1:N) )
!!$    allocate(x_tld_k_of_xk1_dual    (1:N,1:N) )
!!$    allocate(x_tld_k_of_yk1_dual    (1:N,1:N) )
!!$    allocate(x_tld_k_of_zk1_dual    (1:N,1:N) )
!!$    allocate(xp_tld_k_stab_gateaux  (1:N,1:N) )
!!$    allocate(zp_tld_k_stab_gateaux  (1:N,1:N) )
!!$    allocate(xp_tld_k_dual_gateaux  (1:N,1:N) )
!!$    allocate(yp_tld_k_dual_gateaux  (1:N,1:N) )
!!$    allocate(zp_tld_k_dual_gateaux  (1:N,1:N) )
!!$    allocate(FpX_stab               (1:N,1:N) )
!!$    allocate(FpZ_stab               (1:N,1:N) )
!!$    allocate(FpX_dual               (1:N,1:N) )
!!$    allocate(FpZ_dual               (1:N,1:N) )
!!$    allocate(FpY_dual               (1:N,1:N) )
!!$
!!$    !-----------------
!!$    dy_stab     = y_tld_k1_dual - y_k1_dual
!!$    dy_hat_stab = dy_stab/SQRT(SUM(dy_stab**2))
!!$!    WRITE(*,*)' dy_dual ',SQRT(SUM(dy_dual**2))
!!$    dz_stab     = z_tld_k1_stab - z_k1_stab
!!$    dz_hat_stab = dz_stab/SQRT(SUM(dz_stab**2))
!!$!    WRITE(*,*)' dz_stab ',SQRT(SUM(dz_stab**2))
!!$    dx_stab     = x_tld_k1_stab - x_k1_stab
!!$    dx_hat_stab = dx_stab/SQRT(SUM(dx_stab**2))
!!$!    WRITE(*,*)' dx_stab ',SQRT(SUM(dx_stab**2)),SQRT(SUM(x_tld_k1_stab**2)),SQRT(SUM(x_k1_stab**2))
!!$
!!$
!!$    dy_dual     = y_tld_k1_dual - y_k1_dual
!!$    dy_hat_dual = dy_dual/SQRT(SUM(dy_dual**2))
!!$!    WRITE(*,*)' dy_dual ',SQRT(SUM(dy_dual**2))
!!$    dz_dual     = z_tld_k1_dual - z_k1_dual
!!$    dz_hat_dual = dz_dual/SQRT(SUM(dz_dual**2))
!!$!    WRITE(*,*)' dz_dual ',SQRT(SUM(dz_dual**2))
!!$    dx_dual     = x_tld_k1_dual - x_k1_dual
!!$    dx_hat_dual = dx_dual/SQRT(SUM(dx_dual**2))
!!$!    WRITE(*,*)' dx_dual ',SQRT(SUM(dx_dual**2)),SQRT(SUM(x_tld_k1_dual**2)),SQRT(SUM(x_k1_dual**2))
!!$   
!!$    !
!!$    dz_k_dz_k1_stab = MATMUL( dz_hat_stab, m_x_tld_k1_stab )
!!$    dy_k_dy_k1_dual = MATMUL( m_x_tld_k1_dual, dy_hat_dual )
!!$
!!$!    x_tld_k_of_xk1_dual     = MATMUL( MATMUL( m_x_k1_dual,     y_tld_k1_dual ) , MATMUL( z_tld_k1_dual , m_x_k1_dual     )     
!!$    ! dx~_k/dz_k-1 = ( < y~_k | z~_k > - < y~_k | z_k-1 >  m[ x~_k-1 ] )/dz_k-1
!!$    xp_tld_k_stab_gateaux = (x_tld_k_stab &
!!$         -MATMUL(MATMUL(TRANSPOSE(MATMUL(z_tld_k1_stab,m_x_k1_stab )),S_d),MATMUL(z_tld_k1_stab,m_x_k1_stab)) &
!!$                            )/SQRT(SUM(dx_stab**2))
!!$    ! x~'_k_dx_k-1 = m[ x~_k-1 ]^t < dz^t_k-1 | s | z_k >  +  < z_k | s | dz~_k-1 > m[ x~_k-1 ]    
!!$    dz_k_dx_k1_stab = MATMUL(z_tld_k1_stab,mp_stab) 
!!$    FpX_stab        = MATMUL( TRANSPOSE(dz_k_dx_k1_stab) , MATMUL( s_d, z_tld_k_stab ) )  &
!!$                    + MATMUL( TRANSPOSE(z_tld_k_stab)    , MATMUL( s_d, dz_k_dx_k1_stab ) )
!!$    !
!!$    ! dx~_k/dz_k-1 = 
!!$    x_tld_k_of_zk1_stab   = MATMUL( MATMUL( TRANSPOSE(MATMUL(z_k1_stab , m_x_tld_k1_stab )), S_d), MATMUL(z_k1_stab,m_x_tld_k1_stab) )
!!$    zp_tld_k_stab_gateaux = (x_tld_k_stab-x_tld_k_of_zk1_stab)/SQRT(SUM(dz_stab**2))
!!$    ! x~'_k_dz_k-1 = 
!!$    FpZ_stab  = MATMUL(TRANSPOSE(dz_k_dz_k1_stab),MATMUL(s_d,z_tld_k_stab))+MATMUL(TRANSPOSE(z_tld_k_stab),MATMUL(s_d,dz_k_dz_k1_stab))
!!$    !
!!$    !

!!$
!!$    !
!!$    ! dx~_k/dy_k-1 = ( < y~_k | z~_k > - m[ x~_k-1 ] < y_k-1 | z~_k >  )/dy_k-1
!!$    x_tld_k_of_yk1_dual=MATMUL(MATMUL(m_x_tld_k1_dual,y_k1_dual),MATMUL(z_tld_k1_dual,m_x_tld_k1_dual))    
!!$    yp_tld_k_dual_gateaux=(x_tld_k_dual-x_tld_k_of_yk1_dual)/SQRT(SUM(dy_dual**2))
!!$    ! x~'_k_y_k-1 = m[ x~_k-1 ] < dy~_k-1 |z~_k>
!!$    FpY_dual=MATMUL(MATMUL(m_x_tld_k1_dual,dy_hat_dual),z_tld_k_dual) 
!!$    !
!!$    ! dx~_k/dz_k-1 = ( < y~_k | z~_k > - < y~_k | z_k-1 >  m[ x~_k-1 ] )/dz_k-1
!!$    zp_tld_k_dual_gateaux=(x_tld_k_dual-MATMUL( y_tld_k_dual,MATMUL(z_k1_dual,m_x_tld_k1_dual)))/SQRT(SUM(dz_dual**2))
!!$    ! x~'_k_dz_k-1 =  < y~_k | dz^_k-1 > m[ x~_k-1 ] 
!!$    FpZ_dual=MATMUL(y_tld_k_dual,MATMUL(dz_hat_dual,m_x_tld_k1_dual))
!!$    !
!!$    ! dx~_k/dx_k-1 = ( < y~_k | z~_k >  - { m[x_k-1]<y~_k-1|z~_k> + <y~_k|z~_k-1>m[x_k-1] )/dx_k-1
!!$    xp_tld_k_dual_gateaux=(x_tld_k_dual- 0.5d0*( MATMUL(MATMUL(m_x_k1_dual,y_tld_k1_dual),z_tld_k_dual) &
!!$                                               + MATMUL(y_tld_k_dual,MATMUL(z_tld_k1_dual,m_x_k1_dual)) &
!!$                                                ))/SQRT(SUM(dx_dual**2))
!!$    ! x~'_k_dx_k-1 = m[ x~_k-1 ] dx^_k-1 <y~_k-1 | dz~_k > + < y~_k | dz~_k-1 > m[ x~_k-1 ] dx^_k-1
!!$    FpX_dual  = MATMUL( MATMUL( mp_dual , y_tld_k1_dual ), z_tld_k_dual    )  & 
!!$              + MATMUL( y_tld_k_dual , MATMUL( z_tld_k1_dual , mp_dual) )
!!$


!!$
!!$
!!$
!!$!    fpx_dual=(0.5d0)*(dx_hat_dual - MATMUL(x_tld_k_dual,MATMUL(dx_hat_dual,x_tld_k_dual)))
!!$
!!$    FillN_stab=0d0
!!$    FillN_dual=0d0
!!$    do j=1,size(x_tld_k_stab,1)
!!$       FillN_stab=FillN_stab+x_tld_k_stab(j,j)
!!$       FillN_dual=FillN_dual+x_tld_k_dual(j,j)
!!$    enddo
!!$    FillN_stab=(dble(n)-FillN_stab)/dble(n)
!!$    FillN_dual=(dble(n)-FillN_dual)/dble(n)
!!$
!!$    WRITE(*,*)'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
!!$    WRITE(*,*)'   stab: '
!!$    WRITE(*,21)SQRT(SUM(fpx_stab**2)), SQRT(SUM(dx_stab**2))
!!$    WRITE(*,22)SQRT(SUM(fpz_stab**2)), SQRT(SUM(dz_stab**2))
!!$    WRITE(*,24)FillN_stab,0.5d0*(SQRT(SUM(fpx_stab**2))*SQRT(SUM(dx_stab**2))+SQRT(SUM(fpz_stab**2))*SQRT(SUM(dz_stab**2)))
!!$    WRITE(*,*)'   dual: '
!!$    WRITE(*,21)SQRT(SUM(fpx_dual**2)), SQRT(SUM(dx_dual**2))
!!$    WRITE(*,22)SQRT(SUM(fpz_dual**2)), SQRT(SUM(dz_dual**2))
!!$    WRITE(*,23)SQRT(SUM(fpy_dual**2)), SQRT(SUM(dy_dual**2))
!!$    WRITE(*,24)FillN_dual
!!$
!!$21  FORMAT("    f'_dx_k-1 = ",e12.6,", dx_k-1 = ",e12.6)
!!$22  FORMAT("    f'_dz_k-1 = ",e12.6,", dz_k-1 = ",e12.6)
!!$23  FORMAT("    f'_dy_k=1 = ",e12.6,", dy_k-1 = ",e12.6)
!!$24  FORMAT("    [n-trx]/n = ",e12.6,", dx_k   ~ ",e12.6)
!!$
!!$
!!$    WRITE(99,44)kount, tau, FillN_stab, sqrt(sum(dx_stab**2)),sqrt(sum(dz_stab**2)),SQRT(SUM(fpz_stab**2)),SQRT(SUM(fpx_stab**2)), &
!!$                            FillN_dual, sqrt(sum(dx_dual**2)),sqrt(sum(dz_dual**2)),sqrt(sum(dy_dual**2)),SQRT(SUM(fpy_dual**2)), &
!!$                                        SQRT(SUM(fpz_dual**2)),SQRT(SUM(fpx_dual**2))
!!$
!!$44  format(I3,"   ", 20(e12.6,",   "))
!!$
!!$    deallocate(dz_stab)
!!$    deallocate(dz_hat_stab)
!!$    deallocate(dx_stab)
!!$    deallocate(dx_hat_stab)
!!$    deallocate(dy_dual)     
!!$    
!!$    deallocate(dy_hat_dual) 
!!$    deallocate(dz_dual)     
!!$    deallocate(dz_hat_dual) 
!!$    deallocate(dx_dual)     
!!$    deallocate(dx_hat_dual) 
!!$    deallocate(mp_stab) 
!!$    deallocate(mp_dual) 
!!$    deallocate(dz_k_dx_k1_stab)
!!$    deallocate(dz_k_dz_k1_dual)
!!$    deallocate(dy_k_dy_k1_dual)
!!$    deallocate(dz_k_dx_k1_dual)
!!$    deallocate(dy_k_dx_k1_dual)
!!$    deallocate(x_tld_k_of_xk1_stab)
!!$    deallocate(x_tld_k_of_zk1_stab)      
!!$    deallocate(x_tld_k_of_xk1_dual)
!!$    deallocate(x_tld_k_of_yk1_dual)
!!$    deallocate(x_tld_k_of_zk1_dual)   
!!$    deallocate(xp_tld_k_stab_gateaux) 
!!$    deallocate(zp_tld_k_stab_gateaux) 
!!$    deallocate(xp_tld_k_dual_gateaux) 
!!$    deallocate(yp_tld_k_dual_gateaux) 
!!$    deallocate(zp_tld_k_dual_gateaux) 
!!$    deallocate(FpX_stab)  
!!$    deallocate(FpZ_stab)  
!!$    deallocate(FpX_dual)  
!!$    deallocate(FpZ_dual)  
!!$    deallocate(FpY_dual)  
!!$
!!$  END SUBROUTINE SpAMMsand_Error_Analysis2
