module SpAMMsand_inverse_squareroot

  USE spammpack 

  implicit none

  ! Convergence parameters
  REAL(SpAMM_KIND), PARAMETER ::  Approx3  = 2.85d00
  REAL(SPAMM_KIND), PARAMETER ::  ShiftSw  = 5.d-1


contains

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot(s, x, z, tau , t)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: s
    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: x, z, t
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau

    INTEGER                                           :: i, n
    REAL(SpAMM_KIND)                                  :: sc, TrX
    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev

    n=x%frill%ndimn(1)
    FillN=1d10
    !
    z => SpAMM_scalar_plus_tree_2d_symm(SpAMM_one, z)

    DO i = 1, 20

       ! |X_n> = <Z_n|S> |Z_n>
       WRITE(*,99)i,i, i,i, i,i, i,i, i,i, i,i, i,i, i,i, i,i, i,i, i,i
99 format(20(i3))

       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, s, tau*1d-3, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = t )

       write(*,44)' T ',sqrt(t%frill%norm2) ,sqrt(t%child_00%frill%norm2), sqrt(t%child_01%frill%norm2), &
                                             sqrt(t%child_10%frill%norm2), sqrt(t%child_11%frill%norm2)
44     format(A3,'= ',F14.6,', 00=',F14.6,', 01=',F14.6,', 10=',F14.6,', 11=',F14.6)

       x => SpAMM_tree_2d_symm_times_tree_2d_symm( t, z, tau     , alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O = x )

       ! monitor the trace for convergence, maybe look at rate of change at some point too?:
       FillN_prev=FillN
       TrX=SpAMM_trace_tree_2d_symm_recur(x)
       FillN = abs( dble(n) - TrX )/dble(n)       

       WRITE(*,33)tau, i, TrX, FillN
33     format('  ... Tr< ',e6.1,', i=',i2,' > = ', F14.6,' dN=',e10.3)
       !        
       IF(FillN>0.4d0)THEN
          delta=1d-1  ! maybe this should be a variable too, passed in?
          x => spammsand_shift_tree_2d( x, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=spammsand_scaling_invsqrt(SpAMM_zero)
       ELSE
          sc=1d0
       ENDIF

       write(*,*)' x1 = ',sqrt(x%frill%norm2)

       x => spammsand_scaled_invsqrt_mapping( x, sc )

       write(*,*)' x2 = ',sqrt(x%frill%norm2)

       ! |Z_n+1> =  <Z_n| X_n>  
       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=t)

       ! unfortunately, the update could not be done in place 
       z => SpAMM_tree_2d_symm_copy_tree_2d_symm(t, in_O=z)

       write(*,*)' z = ',sqrt(z%frill%norm2)


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

    SHFT=low_new-low_prev*(high_new-low_new)/(high_prev-low_prev)
    SCAL=(high_new-low_new)/(high_prev-low_prev)

    write(*,*)' remapping/shifting  shft = ',shft, ' scal = ',scal


    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    WRITE(*,*)' scal*x = ',SQRT(d%frill%norm2)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)
    WRITE(*,*)' shft(x)= ',SQRT(d%frill%norm2)

  END FUNCTION spammsand_shift_tree_2d

  FUNCTION spammsand_scaled_invsqrt_mapping( x, sc ) result(d) 
    TYPE(spamm_tree_2d_symm), POINTER  :: d
    TYPE(spamm_tree_2d_symm), POINTER  :: x
    REAL(SpAMM_KIND),      INTENT(IN)  :: sc 
    REAL(SpAMM_KIND)                   :: SHFT,SCAL 

    SHFT=SpAMM_half*SQRT(sc)*SpAMM_three
    SCAL=SpAMM_half*(-sc)*SQRT(sc)

    write(*,*)' NS scaling shft = ',shft, ' scal = ',scal

    d => x
    d => SpAMM_scalar_times_tree_2d_symm( scal, d)
    WRITE(*,*)' scal*x = ',SQRT(d%frill%norm2)
    d => SpAMM_scalar_plus_tree_2d_symm(  shft, d)
    WRITE(*,*)' shft(x)= ',SQRT(d%frill%norm2)

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
 
  integer :: i

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)

  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 

  ! the max eigenvalue
  x_hi = SpAMMSand_rqi_extremal(s,1d-8,high_O=.TRUE.)

  ! normalize the max ev of s to 1.  
  s => SpAMM_scalar_times_tree_2d_symm(SpAMM_one/x_hi, s)

  logtau_strt=-8                                       ! starting accuracy
  logtau_stop=-12                                      ! stoping  "
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
