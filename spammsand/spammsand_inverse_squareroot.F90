module SpAMMsand_inverse_squareroot

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

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, t, tau, first)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: x, z, t
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau

    LOGICAL, INTENT(IN) :: first
    INTEGER                                           :: i, n, j
    REAL(SpAMM_KIND)                                  :: sc, TrX, tau_xtra
    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev, &
                                                           trxd, trsd, trtd, trzd

    n=x%frill%ndimn(1)
    FillN=1d10
    !
    z => SpAMM_scalar_plus_tree_2d_symm(SpAMM_one, z)
    
!!$    zd=spamm_zero
!!$    do I=1,s%frill%ndimn(1)
!!$       zd(i,i)=1d0
!!$    enddo
!!$
!!$
!!$    trxd=0d0
!!$    do j=1,n
!!$       trxd=trxd+xd(j,j)
!!$    enddo
!!$
!!$    trsd=0d0
!!$    do j=1,n
!!$       trsd=trsd+sd(j,j)
!!$    enddo
!!$
!!$    WRITE(*,*)' x = ',trxd
!!$    WRITE(*,*)' s = ',trsd

!!$
!!$    CALL SpAMM_convert_tree_2d_symm_to_dense(s,     evec)
!!$    call dsyevd("V", "U", N, evec, N, eval, work, LWORK, iwork, LIWORK, info)
!!$    WRITE(*,*)' EV = ',eval(1),eval(N)


    WRITE(*,*)' s%norm = ',s%frill%norm2

    DO i = 1, 22

       ! |X_n> = <Z_n|S> |Z_n>
!!$       WRITE(*,99)i,i, i,i, i,i, i,i, i,i, i,i, i,i
!!$99     format(20(i3))

!       t => SpAMM_tree_2d_symm_T_times_tree_2d_symm( z, s, tau*1d-3, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = t )

!if(i<3)then

       tau_xtra=tau*1d-1
       if(first)tau_xtra=tau*1d-2 ! stabilize xtra the first step ...

       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s, tau_xtra , NT_O=.FALSE.,alpha_O=SpAMM_zero,beta_O=SpAMM_one, in_O = t )

!else
!       t => SpAMM_tree_2d_symm_times_tree_2d_symm(z,s,tau,NT_O=.FALSE.,alpha_O=SpAMM_zero,beta_O=SpAMM_one, in_O = t )
!endif

!!$!       Td=MATMUL(Zd,Sd) This vs that is a blow up for high cond
!!$       Td=MATMUL(TRANSPOSE(Zd),Sd)
!!$       trtd=0d0
!!$       do j=1,n
!!$          trtd=trtd+td(j,j)
!!$       enddo       
!!$       WRITE(*,*)i,' t = ',trtd       
!!$!       WRITE(*,*)' t  = ',ABS(SQRT(t%frill%norm2)-sqrt(sum(td**2)))/sqrt(sum(td**2))       

       x => SpAMM_tree_2d_symm_times_tree_2d_symm( t, z, tau     , alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = x )

!!$       CALL SpAMM_convert_tree_2d_symm_to_dense(x,     evec)
!!$       call dsyevd("V", "U", N, evec, N, eval, work, LWORK, iwork, LIWORK, info)
!!$       WRITE(*,*)' befor EV = ',eval(1),eval(N)


!!$       Xd=MATMUL(Td,Zd)
!!$       trxd=0d0
!!$       do j=1,n
!!$          trxd=trxd+xd(j,j)
!!$       enddo       
!!$       WRITE(*,*)i,' x = ',trxd      
!!$!       WRITE(*,*)' x  = ',ABS(SQRT(x%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))       
!!$44     format(A3,'= ',F16.8,' dense= ',F16.8)
!!$!       x => SpAMM_convert_dense_to_tree_2d_symm(xd)

       ! monitor the trace for convergence, maybe look at rate of change at some point too?:
       FillN_prev=FillN

       TrX=SpAMM_trace_tree_2d_symm_recur(x)

       if(trx>dble(N).and.first)return

!!$       trxd=0d0
!!$       do j=1,n
!!$       trxd=trxd+xd(j,j)
!!$       enddo

       FillN = abs( dble(n) - TrX )/dble(n)       

!!$       WRITE(*,33)tau, i, TrX, Trxd, FillN
!!$33     format('  ... Tr< ',e6.1,', i=',i2,' > = ', F22.10,F22.10,' dN=',e10.3)

       WRITE(*,33)tau, i, TrX, FillN, z%frill%non0s/dble(N**2)
33     format('  ... Tr< ',e6.1,', i=',i2,' > = ', F18.10,' dN=',e10.3,' Full = ',e5.1)


       !        
       IF(FillN>0.4d0)THEN
          delta=1d-1  ! maybe this should be a variable too, passed in?
!          delta=1d-2  ! maybe this should be a variable too, passed in?
          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_zero)
       ELSE
          sc=1d0
       ENDIF

       x => spammsand_scaled_invsqrt_mapping( x, sc )

!!$       CALL SpAMM_convert_tree_2d_symm_to_dense(x,     evec)
!!$       call dsyevd("V", "U", N, evec, N, eval, work, LWORK, iwork, LIWORK, info)
!!$       WRITE(*,*)' after EV = ',eval(1),eval(N)

!!$       trxd=0d0
!!$       do j=1,n
!!$          trxd=trxd+xd(j,j)
!!$       enddo       
!!$       WRITE(*,*)i,' x = ',trxd

       ! |Z_n+1> =  <Z_n| X_n>  
       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=t)

!!$       td=MATMUL(Zd,Xd)
!!$!       WRITE(*,*)' td  = ',ABS(SQRT(t%frill%norm2)-sqrt(sum(td**2)))/sqrt(sum(td**2))       

       
       ! unfortunately, the update could not be done in place 
       z => SpAMM_tree_2d_symm_copy_tree_2d_symm(t, in_O = z, threshold_O= SQRT(tau) )

!!$       zd=td
!!$       trzd=0d0
!!$       do j=1,n
!!$          trzd=trzd+zd(j,j)
!!$       enddo       
!!$       WRITE(*,*)i,' z = ',trzd
!!$!       WRITE(*,*)' zd  = ',ABS(SQRT(z%frill%norm2)-sqrt(sum(zd**2)))/sqrt(sum(zd**2))       

       IF(FillN<0.1d0.AND.FillN>FillN_prev)THEN
          RETURN
       ELSEIF(FillN<Tau)THEN
          RETURN
       END IF

       ! best acceleration we can hope for
       xo_analytic=xo_analytic*(9d0/4d0)*sc
       !
    END DO
   
  END SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot


  FUNCTION spammsand_shift_tree_2d( x, low_prev, high_prev, low_new, high_new ) RESULT(d)
    !!!!!    shft=low_new+(x-low_prev)*(high_new-low_new)/(high_prev-low_prev)

    TYPE(spamm_tree_2d_symm) ,   POINTER     :: d
    TYPE(spamm_tree_2d_symm) ,   POINTER     :: x
    REAL(SpAMM_KIND), OPTIONAL,INTENT(IN)    :: low_prev, high_prev, low_new, high_new   
    REAL(SpAMM_KIND)                         :: SHFT,SCAL 

    integer :: i

    SHFT=low_new-low_prev*(high_new-low_new)/(high_prev-low_prev)
    SCAL=(high_new-low_new)/(high_prev-low_prev)
!!$!    write(*,*)' remapping/shifting  shft = ',shft, ' scal = ',scal
!!$
    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
!!$    xd=scal*xd
!!$!    WRITE(*,*)' scal*x = ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)
!!$    do i=1,x%frill%ndimn(1)
!!$       xd(i,i)=xd(i,i)+shft
!!$    enddo
!!$    WRITE(*,*)' shft(x)= ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))
  END FUNCTION spammsand_shift_tree_2d

  FUNCTION spammsand_scaled_invsqrt_mapping( x, sc ) result(d) 
    TYPE(spamm_tree_2d_symm), POINTER  :: d
    TYPE(spamm_tree_2d_symm), POINTER  :: x
    REAL(SpAMM_KIND),      INTENT(IN)  :: sc 
    REAL(SpAMM_KIND)                   :: SHFT,SCAL 

    integer :: i

    SHFT=SpAMM_half*SQRT(sc)*SpAMM_three
    SCAL=SpAMM_half*(-sc)*SQRT(sc)
!!$!    write(*,*)' NS scaling shft = ',shft, ' scal = ',scal
    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
!!$    xd=scal*xd
!!$!    WRITE(*,*)' scal*x = ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)
!!$    do i=1,x%frill%ndimn(1)
!!$       xd(i,i)=xd(i,i)+shft
!!$    enddo
!!$!    WRITE(*,*)' shft(x)= ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))
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

  TYPE(spammsand_tree_2d_slices), pointer        :: z, sndwch

  type(SpAMM_tree_2d_symm),       pointer        :: s => null()
  type(SpAMM_tree_2d_symm),       pointer        :: x => null()
  type(SpAMM_tree_2d_symm),       pointer        :: t => null()

!  real(spamm_kind), dimension(:, :), allocatable :: S_dense

  character(len = 1000)                          :: matrix_filename
  real(SpAMM_KIND)                               :: x_hi, x_new, logtau_strt, logtau_stop, logtau_dlta, & 
                                                    tau_dlta, tau_xtra, error, tmp1,tmp2
  logical :: first

  integer, parameter                             :: slices=4

  real(SpAMM_KIND), dimension(1:slices)          :: tau
 
  integer :: i,n,j

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)

  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 

  !=============================================================
  allocate(Sd( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)))
  allocate(Xd( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)))
  allocate(Td( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)))
  allocate(Zd( 1:s%frill%ndimn(1), 1:s%frill%ndimn(2)))

  N = s%frill%ndimn(1)
  LWORK = 1+6*N+2*N**2
  LIWORK = 3+5*N    
  allocate(eval(N))
  allocate(work(LWORK))
  allocate(iwork(LIWORK))

  !  call dsyevd("V", "U", N, X_dense, N, eval, work, LWORK, iwork, LIWORK, info)
  !=============================================================
  ! the max eigenvalue

  x_hi = SpAMMSand_rqi_extremal(s,1d-6,high_O=.TRUE.)

  ! normalize the max ev of s to 1.  
  s => SpAMM_scalar_times_tree_2d_symm(SpAMM_one/x_hi, s)

!!  Sd=s_dense/x_hi
!!  xd=sd

  logtau_strt=-3                                       ! starting accuracy
  logtau_stop=-10                                      ! stoping  "
  logtau_dlta=(logtau_stop-logtau_strt)/dble(slices-1) ! span (breadth) of SpAMM thresholds 
  tau_dlta=10d0**logtau_dlta
    
  allocate(z) 
  sndwch => z ! head of the slices
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
  
  first=.TRUE.
  z=>sndwch

  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>
     
     call spammsand_scaled_newton_shulz_inverse_squareroot( s, x, z%mtx, t, z%tau, first )
     first=.FALSE.
     if(.not.associated(z%nxt))exit

     x => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, in_O = x, threshold_O=1d-12 )
     
     tau_xtra=tau_dlta*z%nxt%tau
     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z%mtx,     x, tau_xtra, NT_O=.FALSE. ,alpha_O=SpAMM_zero,beta_O=SpAMM_one,in_O = t )
     x => SpAMM_tree_2d_symm_times_tree_2d_symm(     t, z%mtx, tau_xtra, NT_O=.TRUE.  ,alpha_O=SpAMM_zero,beta_O=SpAMM_one,in_O = x )

!     CALL SpAMM_convert_tree_2d_symm_to_dense(x,     evec)
!     call dsyevd("V", "U", N, evec, N, eval, work, LWORK, iwork, LIWORK, info)
!     WRITE(*,*)' after EV = ',eval(1),eval(N)

     x_new =  SpAMMSand_rqi_extremal( x, tau_xtra , high_O=.TRUE. )
     x     => SpAMM_scalar_times_tree_2d_symm( SpAMM_one / x_new , x )
     x_hi  = x_hi * x_new

!#     t => SpAMM_tree_2d_symm_times_tree_2d_symm( z%mtx,     x, tau_dlta*z%nxt%tau, NT_O=.TRUE. ,alpha_O=SpAMM_zero,beta_O=SpAMM_one,in_O = t )
!#     x => SpAMM_tree_2d_symm_times_tree_2d_symm(     t, z%mtx, tau_dlta*z%nxt%tau, NT_O=.FALSE. ,alpha_O=SpAMM_zero,beta_O=SpAMM_one,in_O = x )

     s => SpAMM_tree_2d_symm_copy_tree_2d_symm( x, in_O = s , threshold_O = SpAMM_normclean )
     z => z%nxt
     
  enddo

!!$  z=>sndwch
!!$  CALL SpAMM_convert_tree_2d_symm_to_dense(z%mtx,     Z1R)
!!$  z1r=z1r/SQRT(x_hi)
!!$
!!$  CALL SpAMM_convert_tree_2d_symm_to_dense(z%nxt%mtx, Z2R)
!!$
!!$  allocate(Z1L(N,N))
!!$  allocate(Z2L(N,N))
!!$  allocate(x1(N,N))
!!$  allocate(x2(N,N))
!!$  allocate(z_spammsand(N,N))
!!$
!!$
!!$  Z1L=TRANSPOSE(Z1R)
!!$  Z2L=TRANSPOSE(Z2R)
!!$  X1=MATMUL(Z1L,MATMUL(S_dense,Z1R))
!!$  X2=MATMUL(Z2L,MATMUL(X1,Z2R))
!!$
!!$  WRITE(*,*)' nested x '
!!$  WRITE(*,*)''
!!$
!!$  DO J=1,6
!!$     WRITE(*,55)J,(X2(J,I),I=1,6)
!!$  ENDDO
!!$  write(*,*)' '
!!$
!!$  Z_SpAMMsand = matmul(z1r, z2R)
!!$
!!$  WRITE(*,*)' Z_spammsand = '
!!$  WRITE(*,*)''
!!$
!!$  DO J=1,6
!!$     WRITE(*,55)J,(z_spammsand(J,I),I=1,6)
!!$  ENDDO
!!$
!!$
!!$
!!$!  Z_SpAMMsand = (TRANSPOSE(Z_SpAMMsand)+Z_SpAMMsand)*SpAMM_half
!!$  X1=MATMUL(TRANSPOSE(Z_spammsand),MATMUL(S_dense,Z_spammsand))
!!$
!!$  WRITE(*,*)' I_spammsand = '
!!$  WRITE(*,*)''
!!$
!!$  DO J=1,6
!!$     WRITE(*,55)J,(x1(J,I),I=1,6)
!!$55 FORMAT(I4,' ',20(E12.6,', '))
!!$  ENDDO
!!$
!!$  write(*,*)' '
!!$
!!$
!!$  write(*,*)'---------------------------------------------------------- '
!!$


  ! de-normalize, back to max_ev (1->max_ev)...  
  z=>sndwch
  z%mtx => SpAMM_scalar_times_tree_2d_symm(SpAMM_one/SQRT(x_hi), z%mtx)

  ! at this point, we are basically done.  the rest is IO/verification of trace, error, timing & etc
  ! timers down ...

  ! x <= |z_spammsand> = |z_1> . |z_2> ... |z_slices>.  Has to be applied as left (T) and right (N) (its not symmetric)
  x => SpAMM_tree_2d_symm_copy_tree_2d_symm( z%mtx, in_O = x )
  z => z%nxt

  do while(associated(z)) ! build the inverse factors |z> = |z_1>.|z_2> ... |z_slices> (right handed)

     t => SpAMM_tree_2d_symm_times_tree_2d_symm( x, z%mtx, 1d-16, nt_O=.TRUE., & 
                   alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = t )     

     x => SpAMM_tree_2d_symm_copy_tree_2d_symm( t , in_O = x )

     z => z%nxt

  enddo

  ! get the original back ...
  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 
  
  ! I = <z_spammsand|s|z_spammsand>
  t => SpAMM_tree_2d_symm_times_tree_2d_symm( s, x, 1d-16, nt_O=.TRUE. , alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = t )   
  s => SpAMM_tree_2d_symm_times_tree_2d_symm( x, t, 1d-16, nt_O=.FALSE., alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = s )   

  write(*,*)' e1 ',SQRT(ABS(DBLE(N)-s%frill%norm2))/DBLE(N)**2

  ! |error|_F <= [I-1]/N**2
  s => SpAMM_scalar_plus_tree_2d_symm( -SpAMM_one, s) 
  error=SQRT(s%frill%norm2)/dble(N)**2

  write(*,*)' error = ',error

!!$
!!$  evec=S_dense
!!$  call dsyevd("V", "U", N, evec, N, eval, work, LWORK, iwork, LIWORK, info)
!!$
!!$  DO i=1,N
!!$     S_dense(:,I)=evec(:,i)/SQRT(eval(i))
!!$  ENDDO`
!!$
!!$  InvHalf=MATMUL(S_dense,TRANSPOSE(evec))
!!$
!!$  WRITE(*,*)' '
!!$  WRITE(*,*)' invhalf'
!!$  WRITE(*,*)' '
!!$
!!$  DO J=1,6
!!$     WRITE(*,55)J,(invhalf(J,I),I=1,6)
!!$  ENDDO
!!$
!!$
!!$  Error=SUM(ABS(Z_spammsand-InvHalf))/SQRT(SUM(Z_dense**2))
!!$  
!!$  WRITE(*,*)' Error = ',Error

end program SpAMM_sandwich_inverse_squareroot
