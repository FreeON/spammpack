module SpAMMsand_inverse_squareroot

  USE spammpack 

  implicit none

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1

  real(spamm_kind), dimension(:, :), allocatable :: Xd,Zd,Sd,Td
  
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

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, tau , t)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: x, z, t
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau

    INTEGER                                           :: i, n, j
    REAL(SpAMM_KIND)                                  :: sc, TrX
    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev

    n=x%frill%ndimn(1)
    FillN=1d10
    !
    z => SpAMM_scalar_plus_tree_2d_symm(SpAMM_one, z)
    
    zd=spamm_zero
    do I=1,s%frill%ndimn(1)
       zd(i,i)=1d0
    enddo

    WRITE(*,*)'Z0',sqrt(sum(Zd**2))-sum(Zd**2)

    DO i = 1, 20

       ! |X_n> = <Z_n|S> |Z_n>
       WRITE(*,99)i,i, i,i, i,i, i,i, i,i, i,i, i,i
99 format(20(i3))

       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s, tau*1d-3, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = t )

       Td=MATMUL(TRANSPOSE(Zd),Sd)

       WRITE(*,*)' t  = ',ABS(SQRT(t%frill%norm2)-sqrt(sum(td**2)))/sqrt(sum(td**2))       

       x => SpAMM_tree_2d_symm_times_tree_2d_symm( t, z, tau     , alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = x )

       Xd=MATMUL(Td,Zd)
       WRITE(*,*)' x  = ',ABS(SQRT(x%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))       

44     format(A3,'= ',F16.8,' dense= ',F16.8)

!       x => SpAMM_convert_dense_to_tree_2d_symm(xd)

       ! monitor the trace for convergence, maybe look at rate of change at some point too?:
       FillN_prev=FillN


       WRITE(*,*)' FIRST TRACE DONE '

       TrX=SpAMM_trace_tree_2d_symm_recur(x)

       write(*,333)(xd(j,j),j=1,6)
333     format(' 00 = ',10(F12.6,', '))
       write(*,334)(xd(j,j),j=513,519)
334     format(' 11 = ',10(F12.6,', '))



       FillN = abs( dble(n) - TrX )/dble(n)       

       WRITE(*,*)' FIRST TRACE DONE '

       WRITE(*,33)tau, i, TrX, FillN



       trx=0d0
       do j=1,512
          trx=trx+xd(i,i)
          write(*,*)i,trx,xd(i,i)
       enddo
       write(*,*)' tr 00 = ',trx

       trx=0d0
       do j=513,876
          trx=trx+xd(i,i)
          write(*,*)i,trx,xd(i,i)
       enddo
       write(*,*)' tr 11 = ',trx


       WRITE(*,33)tau, i, TrX, FillN
33     format('  ... Tr< ',e6.1,', i=',i2,' > = ', F14.6,' dN=',e10.3)
       !        
       IF(FillN>0.4d0)THEN
!          delta=1d-3  ! maybe this should be a variable too, passed in?
!          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_zero)
       ELSE
          sc=1d0
       ENDIF
       x => spammsand_scaled_invsqrt_mapping( x, sc )

       ! |Z_n+1> =  <Z_n| X_n>  
       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=t)

       td=MATMUL(Zd,Xd)
       
       ! unfortunately, the update could not be done in place 
       z => SpAMM_tree_2d_symm_copy_tree_2d_symm(t, in_O=z)

       zd=td

       WRITE(*,*)' zd  = ',ABS(SQRT(z%frill%norm2)-sqrt(sum(zd**2)))/sqrt(sum(zd**2))       

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

    write(*,*)' remapping/shifting  shft = ',shft, ' scal = ',scal


    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)

    xd=scal*xd
    WRITE(*,*)' scal*x = ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))


    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

    do i=1,x%frill%ndimn(1)
       xd(i,i)=xd(i,i)+shft
    enddo

    WRITE(*,*)' shft(x)= ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))

  END FUNCTION spammsand_shift_tree_2d

  FUNCTION spammsand_scaled_invsqrt_mapping( x, sc ) result(d) 
    TYPE(spamm_tree_2d_symm), POINTER  :: d
    TYPE(spamm_tree_2d_symm), POINTER  :: x
    REAL(SpAMM_KIND),      INTENT(IN)  :: sc 
    REAL(SpAMM_KIND)                   :: SHFT,SCAL 

    integer :: i

    SHFT=SpAMM_half*SQRT(sc)*SpAMM_three
    SCAL=SpAMM_half*(-sc)*SQRT(sc)

    write(*,*)' NS scaling shft = ',shft, ' scal = ',scal

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)

    xd=scal*xd
    WRITE(*,*)' scal*x = ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)

    do i=1,x%frill%ndimn(1)
       xd(i,i)=xd(i,i)+shft
    enddo

    WRITE(*,*)' shft(x)= ',ABS(SQRT(d%frill%norm2)-sqrt(sum(xd**2)))/sqrt(sum(xd**2))
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

  type(SpAMM_tree_2d_symm),       pointer        :: s     => null()
  type(SpAMM_tree_2d_symm),       pointer        :: x_prj => null()
  type(SpAMM_tree_2d_symm),       pointer        :: x_tmp => null()
  real(spamm_kind), dimension(:, :), allocatable :: S_dense
  character(len = 1000)                          :: matrix_filename
  real(SpAMM_KIND)                               :: x_hi, logtau_strt, logtau_stop, logtau_dlta

  integer, parameter                             :: slices=4

  real(SpAMM_KIND), dimension(1:slices)          :: tau
 
  integer :: i,n

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
  x_hi = SpAMMSand_rqi_extremal(s,1d-8,high_O=.TRUE.)

  ! normalize the max ev of s to 1.  
  s => SpAMM_scalar_times_tree_2d_symm(SpAMM_one/x_hi, s)
  Sd=s_dense/x_hi

  logtau_strt=-16                                       ! starting accuracy
  logtau_stop=-32                                      ! stoping  "
  logtau_dlta=(logtau_stop-logtau_strt)/dble(slices-1) ! span (breadth) of SpAMM thresholds 
    
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
  x_prj => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
  x_tmp => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
  
  z=>sndwch
  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>
     
     call spammsand_scaled_newton_shulz_inverse_squareroot( s, x_prj, z%mtx, z%tau, x_tmp )
     
     if(.not.associated(z%nxt))exit

     x_prj => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, x_prj )


     x_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( x_prj,  z%mtx, z%nxt%tau, &
                                                         alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=x_tmp)
     
     x_prj => SpAMM_tree_2d_symm_times_tree_2d_symm( z%mtx,  x_tmp, z%nxt%tau, &
                                                         alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=x_prj)     

     s => SpAMM_tree_2d_symm_copy_tree_2d_symm( x_prj, s )

     z=>z%nxt

     write(*,*)' '
     
  enddo

  ! de-normalize the 1 to max_ev.  
  sndwch%mtx => SpAMM_scalar_times_tree_2d_symm(x_hi, sndwch%mtx)

  ! check the result 



end program SpAMM_sandwich_inverse_squareroot
