module SpAMMsand_inverse_squareroot
! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan  -Wextra -std=f2008 "
  USE spammpack 

  implicit none

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1

  real(spamm_kind), dimension(:, :), allocatable :: Xd,Zd,Sd,Td
  real(spamm_kind), dimension(:, :), allocatable :: Z_dense, InvHalf, S_dense,z1,z2,z3,z_spammsand


  real(spamm_kind), dimension(:, :), allocatable :: Z1L,Z2L,Z1R,Z2R,X1,X2 

  
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
  
contains

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, t, tau, first, second, kount)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER :: x, z, y, t
!    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: x, z, t
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau
!    LOGICAL, INTENT(IN)                               :: first
    LOGICAL                              :: first,second
    INTEGER, INTENT(INOUT)                            :: kount
    INTEGER                                           :: i, n, j
    REAL(SpAMM_KIND)                                  :: sc, TrX, tau_xtra

    REAL(SpAMM_KIND)                                  :: twist_x,twist_z,twist_y,twist_zs
    REAL(SpAMM_KIND)                                  :: x_work,z_work,y_work,zs_work

    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev, &
                                                         trxd, trsd, trtd, trzd

    real(spamm_kind), dimension(:,:), allocatable :: x_d
    n=x%frill%ndimn(1)
    FillN=1d10
    
!    first=.false.

    IF(.nOt.First)tHeN
       ! y_0 => s
       y => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
       y => SpAMM_tree_2d_symm_copy_tree_2d_symm( S , in_o = y, threshold_O = SpAMM_normclean ) 
    endIF

    ! z_0 => I
    z=>SpAMM_set_identity_2d_symm ( s%frill%ndimn, in_o = z )

    !  x_0 => s
    x => SpAMM_tree_2d_symm_copy_tree_2d_symm( s , in_O = x, threshold_O = SpAMM_normclean ) 

    DO i = 0, 26

       ! accounting
       IF(first.and.i==0)then
          kount=1
       ELSE
          kount=kount+1
       ENDIF

       ! monitor the trace for convergence, maybe look at rate of change at some point too?:
       FillN_prev=FillN
       TrX=SpAMM_trace_tree_2d_symm_recur(x)
       FillN = abs( dble(n) - TrX )/dble(n)       

!       write(*,*)' trx = ',trx, ' n  = ',n,' filln = ',filln

       WRITE(101,*)kount, filln
       WRITE(*,33)tau, kount, TrX, FillN, 1d2*x%frill%non0s/dble(x%frill%ndimn(1))**2, y_work*1d2 , z_work*1d2 , zs_work*1d2 , x_work*1d2 
!twist_y,twist_z,twist_zs,twist_x

       IF(i>2 .and. FillN<0.1d0 .AND. FillN>FillN_prev )then
!          WRITE(*,*)' fill n = ',filln,' filln_prev ',filln_prev
!          write(*,*)' elevation' 
          RETURN  ! Elevation
       endif
       IF( FillN <  Tau*1d1 )then
!          WRITE(*,*)' fill n = ',filln,' tau**2 = ',tau
!          write(*,*)' anhiliaiton '
          RETURN  ! Anihilation
       end IF

33     format('  ... Tr< ',e6.1,', i=',i2,' > = ', F18.10,' dN=',e10.3,' %of N^2 = ',f10.5,'%,  %of N^3 = ',f10.5,'%, ',f10.5,'%, ',f10.5,'%, ',f10.5,'%')
!e10.3,', ',e10.3,', ',e10.3,', ',e10.3)


       IF(FillN>0.4d0)THEN
          delta=1.0d-1 ! maybe this should be a variable too, passed in?
          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_zero)
       ELSEIF(FillN>0.1d0)THEN
          delta=1.0d-2 ! maybe this should be a variable too, passed in?
          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_half)
       ELSE
          sc=1d0
       ENDIF

!!$       !=====================================================================
!!$       CALL SpAMM_convert_tree_2d_symm_to_dense(x,xd)
!!$       xd=SpAMM_half*(xd+TRANSPOSE(xd))
!!$       x => SpAMM_convert_dense_to_tree_2d_symm( xd, in_O = x )
!!$       !=====================================================================

!!$       call dsyevd("V", "U", N, Xd, N, eval, work, LWORK, iwork, LIWORK, info)


!!$       call dsyevd("V", "U", N, Xd, N, eval, work, LWORK, iwork, LIWORK, info)
!!$       sc=spammsand_scaling_invsqrt(eval(1))
!!$       WRITE(*,*)' eval = ',eval(1),eval(n),' sc = ',sc
!!$
!!       sc=1d0
!!$
       x => spammsand_scaled_invsqrt_mapping( x, sc )

       ! |Z_n+1> =  <Z_n| X_n>  
       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau , nt_O=.TRUE., in_O = t )
       z_work=t%frill%flops/dble(t%frill%ndimn(1))**3

       ! update (threshold)
       z => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = z, threshold_O = tau )



       !=====================================================================
!!$       CALL SpAMM_convert_tree_2d_symm_to_dense(z,xd)
!!$       xd=SpAMM_half*(xd+TRANSPOSE(xd))
!!$       z => SpAMM_convert_dense_to_tree_2d_symm( xd, in_O = z )
       !=====================================================================

!       CALL SpAMM_convert_tree_2d_symm_to_dense(z, twist)
!       twist_z=SQRT(SUM(twist-TRANSPOSE(twist))**2)
!       twist_z = SpAMM_twist_tree_2d_symm(z)

       IF(.nOt.First)tHeN
          ! <Y_n+1|
          t => SpAMM_tree_2d_symm_times_tree_2d_symm( x, y, tau , nt_O=.TRUE., in_O = t )
          y_work=t%frill%flops/dble(t%frill%ndimn(1))**3
          ! update (threshold)
          y => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = y, threshold_O = tau ) 

          !=====================================================================
!!$          CALL SpAMM_convert_tree_2d_symm_to_dense(y,xd)
!!$          xd=SpAMM_half*(xd+TRANSPOSE(xd))
!!$          y => SpAMM_convert_dense_to_tree_2d_symm( xd, in_O = y )
          !=====================================================================

!          twist_y =SQRT( SpAMM_twist_2d_symm_recur (y) 
!          CALL SpAMM_convert_tree_2d_symm_to_dense(y, twist)
!          twist_y=SQRT(SUM(twist-TRANSPOSE(twist))**2)

       ENDIF

       IF(first)then

          ! <X_n> = <Z_n|S|Z_n>
          if(second)then
             tau_xtra=tau*1d-2  ! xtra stabilization on first multiply 
          else
             tau_xtra=tau       
          endif

          write(*,*)' tau extra = ',tau_xtra

!          t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s, tau  , NT_O=.false. , in_O = t )
!          t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s, tau  , NT_O=.TRUE. , in_O = t )

          t => SpAMM_tree_2d_symm_times_tree_2d_symm( s, z, tau  , NT_O=.TRUE. , in_O = t )

          WRITE(*,*)'BEFORE S_OT_Z ',t%frill%norm2


          zs_work=t%frill%flops/dble(t%frill%ndimn(1))**3

          !=====================================================================
!!$          CALL SpAMM_convert_tree_2d_symm_to_dense(t,xd)
!!$          xd=SpAMM_half*(xd+TRANSPOSE(xd))
!!$          t => SpAMM_convert_dense_to_tree_2d_symm( xd, in_O = t )
          !=====================================================================

!          write(*,*)' z norm = ',SQRT(z%frill%norm2),' total work = ',t%frill%flops/dble(t%frill%ndimn(1))**3

!          x => SpAMM_tree_2d_symm_times_tree_2d_symm( t, z, tau  , NT_O=.TRUE.  , in_O = x )

!          x => SpAMM_tree_2d_symm_times_tree_2d_symm( z, t,  tau  , NT_O=.FALSE. , in_O = x )

          x => SpAMM_tree_2d_symm_times_tree_2d_symm( z, t,  tau  , NT_O=.TRUE. , in_O = x )

          CALL SpAMM_convert_tree_2d_symm_to_dense( x, x_d )

          WRITE(*,*)' XPERP = ',SQRT(SUM( ( x_d-TRANSPOSE(x_d))**2))

          IF(I==6)CALL SpAMM_demo_spamm_stabilized_factorization(z,s, tau)

          x_work=x%frill%flops/dble(x%frill%ndimn(1))**3

          !=====================================================================
!!$          CALL SpAMM_convert_tree_2d_symm_to_dense(x,xd)
!!$          xd=SpAMM_half*(xd+TRANSPOSE(xd))
!!$          x => SpAMM_convert_dense_to_tree_2d_symm( xd, in_O = x )
          !=====================================================================
          
!          write(*,*)' x norm = ',SQRT(x%frill%norm2),' total work = ',x%frill%flops/dble(x%frill%ndimn(1))**3
!          CALL SpAMM_convert_tree_2d_symm_to_dense(t, twist)
!          twist_zs=SQRT(SUM(twist-TRANSPOSE(twist))**2)
!          CALL SpAMM_convert_tree_2d_symm_to_dense(x, twist)
!          twist_x=SQRT(SUM(twist-TRANSPOSE(twist))**2)
!          twist_zs=SQRT( SpAMM_twist_2d_symm_recur (t) )
!          twist_x =SQRT( SpAMM_twist_2d_symm_recur (x) )

       else
          ! <X_n> = <Y_n|Z_n>
          x => SpAMM_tree_2d_symm_times_tree_2d_symm( y, z, tau   , nt_O=.TRUE. , in_O = x )
          x_work=x%frill%flops/dble(x%frill%ndimn(1))**3

          !=====================================================================
!!$          CALL SpAMM_convert_tree_2d_symm_to_dense(x,xd)
!!$          xd=SpAMM_half*(xd+TRANSPOSE(xd))
!!$          x => SpAMM_convert_dense_to_tree_2d_symm( xd, in_O = x )
          !=====================================================================

!          CALL SpAMM_convert_tree_2d_symm_to_dense(x, twist)
!          twist_x=SQRT(SUM(twist-TRANSPOSE(twist))**2)
!          twist_x =SQRT( SpAMM_twist_2d_symm_recur (x) )

       endif
!!$
!!$       ! best acceleration we can hope for
!!$       xo_analytic=xo_analytic*(9d0/4d0)*sc
       !

    END DO

    IF(.not.first)call SpAMM_destruct_tree_2d_symm_recur (y)
   
  END SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot


  SUBROUTINE SpAMM_demo_spamm_stabilized_factorization(z,s,tau)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)   :: z,s

    TYPE(SpAMM_tree_2d_symm), POINTER :: z_sym=>NULL(),s_ot_z=>NULL(), &
         zT_ot_s_ot_z=>NULL(), z_ot_s_ot_z=>NULL(), zT=>NULL(), s_ot_zT=>NULL()

    real(spamm_kind), dimension(:,:), allocatable    :: &
    s_d, z_d, x_stab, x_naiv, s_dot_z_d,  s_dot_zT_d,  &
    s_ot_z_d, s_ot_zT_d, delta_s_ot_z_d, delta_s_ot_zT_d,  &
    delta_zT_ot_s_ot_z_d,e_perp,delta_z_ot_s_ot_z_d, &
    z_dot_delta_s_ot_z_d,       zT_dot_delta_s_ot_zT_d,zT_dot_delta_s_ot_z_d,zT_dot_deltaT_s_ot_z_d, &
    z_ot_s_ot_z_d,zT_ot_s_ot_z_d

    INTEGER :: M

    REAL(SpAMM_kind), intent(in) ::  tau 

    M=s%frill%ndimn(1)

    ALLOCATE(s_d                  (1:M,1:M))
    ALLOCATE(z_d                  (1:M,1:M))
    ALLOCATE(x_naiv               (1:M,1:M))
    ALLOCATE(x_stab               (1:M,1:M))
    ALLOCATE(e_perp               (1:M,1:M))
    ALLOCATE(            s_ot_z_d (1:M,1:M))
    ALLOCATE(           s_dot_z_d (1:M,1:M))
    ALLOCATE(      delta_s_ot_z_d (1:M,1:M))
    ALLOCATE(delta_z_ot_s_ot_z_d  (1:M,1:M))
    ALLOCATE(delta_zT_ot_s_ot_z_d (1:M,1:M))
    ALLOCATE(       z_ot_s_ot_z_d (1:M,1:M))
    ALLOCATE(z_dot_delta_s_ot_z_d (1:M,1:M))
    ALLOCATE(zT_dot_delta_s_ot_z_d(1:M,1:M))
    ALLOCATE(zT_dot_deltaT_s_ot_z_d(1:M,1:M))

    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d )
    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_d )
    zT => SpAMM_convert_dense_to_tree_2d_symm( TRANSPOSE(z_d), in_O = zT )

    ! Symmetrize 

    s_dot_z_d =MATMUL( s_d , z_d )
    s_dot_zT_d=MATMUL( s_d , TRANSPOSE(z_d) )

    ! Errors are relative to a targeted functional form: <zt|s|z> or <z|s|z>  
    ! for some incoming factor z, which is hopefully the result of stable iteration 
    x_naiv   =MATMUL(           z_d  , s_dot_z_d ) 
    x_stab   =MATMUL( TRANSPOSE(z_d) , s_dot_z_d ) 
 
    ! The approximate (ot is the spamm multiply, dot the matmul) intermediate and final NS products:
          s_ot_z => SpAMM_tree_2d_symm_times_tree_2d_symm( s,      z, tau  , NT_O=.TRUE.  , in_O =       s_ot_z )
         s_ot_zT => SpAMM_tree_2d_symm_times_tree_2d_symm( s,     zT, tau  , NT_O=.TRUE.  , in_O =      s_ot_zT )
     z_ot_s_ot_z => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s_ot_z, tau  , NT_O=.TRUE.  , in_O =  z_ot_s_ot_z )
    zT_ot_s_ot_z => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s_ot_z, tau  , NT_O=.FALSE. , in_O = zT_ot_s_ot_z )
    CALL SpAMM_convert_tree_2d_symm_to_dense(       s_ot_z  ,       s_ot_z_d  )
    CALL SpAMM_convert_tree_2d_symm_to_dense(       s_ot_zT ,       s_ot_zT_d )
    CALL SpAMM_convert_tree_2d_symm_to_dense(  z_ot_s_ot_z  ,  z_ot_s_ot_z_d  )
    CALL SpAMM_convert_tree_2d_symm_to_dense( zT_ot_s_ot_z  , zT_ot_s_ot_z_d  )

    e_perp=z_ot_s_ot_z_d-TRANSPOSE(z_ot_s_ot_z_d)
    write(*,*)' perp norm = ',SQRT(SUM(e_perp**2))

    e_perp=zT_ot_s_ot_z_d-TRANSPOSE(zT_ot_s_ot_z_d)
    write(*,*)' perp norm = ',SQRT(SUM(e_perp**2))

    ! Here are primary errors (first multiply)
            delta_s_ot_z_d  = s_ot_z_d  - s_dot_z_d  
            delta_s_ot_zT_d = s_ot_zT_d - s_dot_zT_d  



       z_dot_delta_s_ot_z_d  = MATMUL(          z_d ,           delta_s_ot_z_d )
      zT_dot_delta_s_ot_zT_d = MATMUL(TRANSPOSE(z_d) ,          delta_s_ot_zT_d )
      zT_dot_delta_s_ot_z_d  = MATMUL(TRANSPOSE(z_d),           delta_s_ot_z_d )
     zT_dot_deltaT_s_ot_z_d  = MATMUL(TRANSPOSE(z_d), TRANSPOSE(delta_s_ot_z_d))

    ! Here are secondary errors
     delta_z_ot_s_ot_z_d =  z_ot_s_ot_z_d - x_stab -  z_dot_delta_s_ot_z_d
    delta_zT_ot_s_ot_z_d = zT_ot_s_ot_z_d - x_stab - zT_dot_delta_s_ot_z_d


    WRITE(*,*)' Check 1 A',SQRT(SUM( (s_ot_z_d   - (s_dot_z_d  + delta_s_ot_z_d)  )**2))
    WRITE(*,*)' Check 1 B ',SQRT(SUM( (s_ot_zT_d - (s_dot_zT_d + delta_s_ot_zT_d) )**2))




    write(*,*)' ||z.delta_s_ot_z||_F    = ',SQRT(SUM( z_dot_delta_s_ot_z_d**2 )), &
         ' < ',Tau*SQRT(s%frill%norm2)*z%frill%norm2

    write(*,*)' naive perp error        = ',SQRT(SUM( (z_dot_delta_s_ot_z_d  &
                                           - TRANSPOSE(z_dot_delta_s_ot_z_d) )**2 ))
!                                                     -zT_dot_delta_s_ot_zT_d  )**2 ))

    write(*,*)' stabilized perp error   = ',SQRT(SUM( (zT_dot_delta_s_ot_z_d  &
                                           - TRANSPOSE(zT_dot_delta_s_ot_z_d) )**2 ))


    write(*,*)' naive perp error        = ',SQRT(SUM( (delta_z_ot_s_ot_z_d   &
                                           - TRANSPOSE(delta_z_ot_s_ot_z_d) )**2 ))

    write(*,*)' stabilized perp error   = ',SQRT(SUM( (delta_zT_ot_s_ot_z_d  &
                                           - TRANSPOSE(delta_zT_ot_s_ot_z_d) )**2 ))
    DEALLOCATE(s_d                  )
    DEALLOCATE(z_d                  )
    DEALLOCATE(x_stab               )
    DEALLOCATE(e_perp               )
    DEALLOCATE(s_ot_z_d             )
    DEALLOCATE(s_dot_z_d            )
    DEALLOCATE(delta_s_ot_z_d       )
    DEALLOCATE(delta_z_ot_s_ot_z_d  )
    DEALLOCATE(delta_zT_ot_s_ot_z_d )
    DEALLOCATE(zT_dot_deltaT_s_ot_z_d)

    stop

  END SUBROUTINE SpAMM_demo_spamm_stabilized_factorization


  SUBROUTINE SpAMM_demo_spamm_stabilized_factorization2(z,s,tau)
    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)   :: z,s
    TYPE(SpAMM_tree_2d_symm) , POINTER               :: z_symm=>NULL(),x_tilde=>NULL(),sz_tilde=>NULL(), &
         zsz_tilde=>NULL(), sz=>NULL(), z_ot_s_dot_z_tilde=>NULL()
    real(spamm_kind), dimension(:,:), allocatable    :: s_d, z_d, x_d,x_tilde_d, e_perp, &
         delta_sz_d, delta_zsz_d,z_dot_delta_sz_d, z_dot_sz_tilde_d, sz_tilde_d,z_ot_s_dot_z_d, &
         s_dot_z_d, sz_d, zsz_tilde_d
    INTEGER :: M

    REAL(SpAMM_kind), intent(in) ::  tau 

    sz_tilde => SpAMM_tree_2d_symm_times_tree_2d_symm( s,         z, tau*1d-2  , NT_O=.TRUE. , in_O = sz_tilde )
    x_tilde  => SpAMM_tree_2d_symm_times_tree_2d_symm( z, sz_tilde , tau  , NT_O=.TRUE. , in_O = x_tilde )

    M=s%frill%ndimn(1)

    ALLOCATE(x_d(1:M,1:M))
    ALLOCATE(s_d(1:M,1:M))
    ALLOCATE(sz_d(1:M,1:M))
    ALLOCATE(x_tilde_d(1:M,1:M))
    ALLOCATE(e_perp(1:M,1:M))
    ALLOCATE(s_dot_z_d(1:M,1:M))
    ALLOCATE(delta_sz_d(1:M,1:M))
    ALLOCATE(delta_zsz_d(1:M,1:M))
    ALLOCATE(z_dot_delta_sz_d(1:M,1:M))
    ALLOCATE(sz_tilde_d(1:M,1:M))
    ALLOCATE(zsz_tilde_d(1:M,1:M))
    ALLOCATE(z_ot_s_dot_z_d(1:M,1:M))

    CALL SpAMM_convert_tree_2d_symm_to_dense( s, s_d )
    CALL SpAMM_convert_tree_2d_symm_to_dense( z, z_d )

    s_dot_z_d=MATMUL( s_d , z_d )
    x_d=MATMUL(z_d,s_dot_z_d)

    write(*,*)' z perp ',SQRT(SUM( ( z_d-TRANSPOSE(z_d) )**2))
 
    sz => SpAMM_convert_dense_to_tree_2d_symm( s_dot_z_d , in_O = sz )
    CALL SpAMM_convert_tree_2d_symm_to_dense( sz_tilde, sz_tilde_d )
    delta_sz_d = sz_tilde_d - sz_d  
    WRITE(*,*)' First check ',SQRT(SUM( (sz_tilde_d - (sz_d+delta_sz_d) )**2))

    ! note tilde is for z \ot [s.z]:
    z_ot_s_dot_z_tilde => SpAMM_tree_2d_symm_times_tree_2d_symm( z, sz , tau  , NT_O=.TRUE. , in_O = zsz_tilde )
    CALL SpAMM_convert_tree_2d_symm_to_dense( zsz_tilde, z_ot_s_dot_z_d )

    z_dot_delta_sz_d = MATMUL( z_d , delta_sz_d )    ! note tilde is for z \ot [s.z]:
    zsz_tilde => SpAMM_tree_2d_symm_times_tree_2d_symm( z, sz , tau  , NT_O=.TRUE. , in_O = zsz_tilde )

    delta_zsz_d = (zsz_tilde_d - x_d) - z_dot_delta_sz_d 
    WRITE(*,*)' Second check ',SQRT(SUM( ( zsz_tilde_d - (x_d+delta_zsz_d +z_dot_delta_sz_d))**2 ))

    e_perp=x_d-TRANSPOSE(x_d)
    write(*,*)' perp error = ',SQRT(SUM(e_perp**2))

    e_perp = delta_zsz_d + z_dot_delta_sz_d - TRANSPOSE( delta_zsz_d + z_dot_delta_sz_d )
    write(*,*)' perp error = ',SQRT(SUM(e_perp**2))
    
    write(*,*)' ||z.delt_sz||_F  ',SQRT(SUM(z_dot_delta_sz_d)**2),' < ',Tau*SQRT(s%frill%norm2)*z%frill%norm2


    write(*,*)' ||delta_zsz||_F = ',SQRT(SUM(delta_zsz_d**2)),'<',Tau*(1d0+tau)*SQRT(s%frill%norm2)*z%frill%norm2



    e_perp=MATMUL( TRANSPOSE(delta_sz_d),z_d)-MATMUL( TRANSPOSE(z_d), delta_sz_d)
    write(*,*)' perp error = ',SQRT(SUM(e_perp**2))


    call SpAMM_destruct_tree_2d_symm_recur (x_tilde)
    call SpAMM_destruct_tree_2d_symm_recur (sz_tilde)

    DEALLOCATE(x_d)
    DEALLOCATE(s_d)
    DEALLOCATE(e_perp)
    DEALLOCATE(x_tilde_d)
    DEALLOCATE(sz_tilde_d)
    DEALLOCATE(zsz_tilde_d)
    DEALLOCATE(delta_sz_d)
    DEALLOCATE(z_dot_delta_sz_d)
    DEALLOCATE(z_ot_s_dot_z_d)
    DEALLOCATE(s_dot_z_d)


  END SUBROUTINE SpAMM_demo_spamm_stabilized_factorization2

!!$



  FUNCTION spammsand_shift_tree_2d( x, low_prev, high_prev, low_new, high_new ) RESULT(d)
    !!!!!    shft=low_new+(x-low_prev)*(high_new-low_new)/(high_prev-low_prev)

    TYPE(spamm_tree_2d_symm) ,   POINTER     :: d
    TYPE(spamm_tree_2d_symm) ,   POINTER     :: x
    REAL(SpAMM_KIND), OPTIONAL,INTENT(IN)    :: low_prev, high_prev, low_new, high_new   
    REAL(SpAMM_KIND)                         :: SHFT,SCAL 

    integer :: i

    SHFT=low_new-low_prev*(high_new-low_new)/(high_prev-low_prev)
    SCAL=(high_new-low_new)/(high_prev-low_prev)

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
 
  integer :: i,n,j,k, kount

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)
  S_dense=SpAMM_half*(S_dense+TRANSPOSE(S_dense))


  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm( S_DENSE, in_O = s )

  !=============================================================
!!$  allocate(Sd( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)) )
!!$  allocate(Xd( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)) )
!!$  allocate(Td( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)) )
!!$  allocate(Zd( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)) )
!!$
!!$  N = s%frill%ndimn(1)
!!$  LWORK = 1+6*N+2*N**2
!!$  LIWORK = 3+5*N    
!!$  allocate(eval(N))
!!$  allocate(work(LWORK))
!!$  allocate(iwork(LIWORK))

  !  call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
  !=============================================================
  ! the max eigenvalue

  x_hi = SpAMMSand_rqi_extremal(s,1d-4,high_O=.TRUE.)
  WRITE(*,*)' hi extremal = ',x_hi

  ! normalize the max ev of s to 1.  
  s       => SpAMM_scalar_times_tree_2d_symm( SpAMM_one/x_hi , s )
  s_orgnl => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, in_O = s_orgnl, threshold_O = SpAMM_normclean )

  logtau_strt=-3                                       ! starting accuracy
  logtau_stop=-9                                        ! stoping  "
  logtau_dlta=(logtau_stop-logtau_strt)/dble(slices-1) ! span (breadth) of SpAMM thresholds 
  tau_dlta=10d0**logtau_dlta
  final_tau=10d0**logtau_stop                          ! penultimate spamm threshold
    
  allocate(z) 
  z_head => z ! head of the slices
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

  kount=0
  first=.TRUE.
  second=.TRUE.

  z=>z_head
  ! z_total => I
  z_total=>SpAMM_set_identity_2d_symm ( s%frill%ndimn, in_o = z_total )

  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>
     
     call spammsand_scaled_newton_shulz_inverse_squareroot( s, x, z%mtx, t, z%tau, first, second, kount )

     first=.FALSE.
     second=.FALSE.
     if(.not.associated(z%nxt))exit    

     tau_xtra=z%tau
!     tau_xtra=1d-3*z%tau

     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z_total, z%mtx, z%nxt%tau , NT_O=.TRUE., in_O = t ) !z%nxt%tau, NT_O=.FALSE., in_O = t )
     z_total => SpAMM_tree_2d_symm_copy_tree_2d_symm( t, in_O = z_total, threshold_O = z%tau )       ! *SQRT(SQRT(z%mtx%frill%norm2)) )

     ! Nested sandwich, to recompute the error every time. 
     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z_total, s_orgnl, z%nxt%tau, NT_O=.FALSE., in_O = t ) !z%nxt%tau, NT_O=.FALSE., in_O = t )
     s => SpAMM_tree_2d_symm_times_tree_2d_symm(       t, z_total, z%nxt%tau, NT_O=.TRUE. , in_O = s )

!     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z%mtx,     s, final_tau, NT_O=.FALSE., in_O = t ) !z%nxt%tau, NT_O=.FALSE., in_O = t )
!     s => SpAMM_tree_2d_symm_times_tree_2d_symm(     t, z%mtx, final_tau, NT_O=.TRUE. , in_O = s )

      s_work=s%frill%flops/dble(s%frill%ndimn(1))**3
     zs_work=t%frill%flops/dble(t%frill%ndimn(1))**3


!!$     x_new =  SpAMMSand_rqi_extremal( s, 1d-8 , high_O=.TRUE. )
!!$     WRITE(*,*)' hi extremal = ',x_new 
!!$
!!$     x     => SpAMM_scalar_times_tree_2d_symm( SpAMM_one / x_new , x )
!!$     x_hi  = x_hi * x_new
!!$
!!$


     write(*,*)' t work = ',zs_work,', s_work = ',s_work

     z => z%nxt


!!$     x_new =  SpAMMSand_rqi_extremal( x, z%nxt%tau , high_O=.TRUE. )
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
