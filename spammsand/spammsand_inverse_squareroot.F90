module SpAMMsand_inverse_squareroot

  USE spammpack 

  implicit none

contains

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot( x, z, tau , t)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: x, z, t
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau

    INTEGER                                           :: i, n
    REAL(SpAMM_KIND)                                  :: sc
    REAL(SpAMM_KIND)                                  :: xo_analytic, delta, FillN, FillN_prev

    FillN=1d10
    !
    DO i = 1, 22

       ! |X_n> = <Z_n|S> |Z_n>
       t => SpAMM_tree_2d_symm_times_tree_2d_symm(z,x,tau*1d-3,alpha_O=SpAMM_zero,beta_O=SpAMM_one,in_O=t)
       x => SpAMM_tree_2d_symm_times_tree_2d_symm(t,z,tau     ,alpha_O=SpAMM_zero,beta_O=SpAMM_one,in_O=x)
       
       ! monitor the trace for convergence, maybe look at rate of change at some point too?:
       FillN_prev=FillN
       FillN = ( dble(n) - SpAMM_trace_tree_2d_symm_recur(x) )/dble(n)       

       !        
       IF(FillN>0.4d0)THEN
          delta=1d-1  ! maybe this should be a variable too, passed in?
!          X => GSOLVE_SPECTRAL_Shift( X, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
!          sc=ScaleInvSqrt(0d0)
       ELSE
!          sc=1d0
       ENDIF

!       X =>  GSOLVE_SPECTRAL_InvSqrt_Scaled_NS( X, sc ) 

       ! |Z_n+1> =  <Z_n| X_n>  
       t => SpAMM_tree_2d_symm_times_tree_2d_symm( z, x, tau, alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=t)
       ! unfortunately, the update could not be done in place 
       z => SpAMM_tree_2d_symm_copy_tree_2d_symm(t, in_O=z)
     
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

end module SpAMMsand_inverse_squareroot



! Nested SpAMM solvers (SpAMM sandwitches) for matrix functions
program SpAMM_sandwich_inverse_squareroot

  USE spammpack 
  USE sandpack 
  USE SpAMMsand_inverse_squareroot
  use test_utilities
  
  implicit none


  TYPE(spammsand_tree_2d_slices), pointer        :: z, sandwtch
  type(SpAMM_tree_2d_symm),       pointer        :: S     => null()
  type(SpAMM_tree_2d_symm),       pointer        :: x_prj => null()
  type(SpAMM_tree_2d_symm),       pointer        :: x_tmp => null()
  real(spamm_kind), dimension(:, :), allocatable :: S_dense
  character(len = 1000)                          :: matrix_filename
  real(SpAMM_KIND)                               :: x_hi, logtau_strt, logtau_stop, logtau_dlta

  integer :: i

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)

  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 

  ! the max eigenvalue
  x_hi=SpAMMSand_rqi_extremal(s,1d-4,high_O=.TRUE.)

  ! normalize the max ev of s to 1.  
  s => SpAMM_scalar_times_tree_2d_symm(SpAMM_one/x_hi, s)

  logtau_strt=-4
  logtau_stop=-12
  logtau_dlta=(logtau_stop-logtau_strt)/dble(slices+1) 
    
  sandwtch => z ! head of the slices

  do i=1,slices

     allocate(z) 

     z%tau =  10d0**( logtau_strt + logtau_dlta*float(i-1) )
     z%mtx => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
     WRITE(*,*)' tau = ',z%tau, sandwtch%tau

     z     => z%nxt

  enddo

  ! the projector and its residual slices
  x_prj => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, x_prj )

  ! temporary work space ...
  x_tmp => SpAMM_new_top_tree_2d_symm( s%frill%ndimn )
  
  z=>sandwtch
  do while(associated(z)) ! build the nested inverse factors |z> = |z_1>.|z_2> ... |z_s>

!     call spammsand_scaled_newton_shulz_inverse_squareroot( x_nxt, z%mtx(i), z%tau(i), x_tmp )

     if(.not.associated(z%nxt))exit

!     x_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( x_prj,  z%mtx, z%nxt%tau, &
!                                                     alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=x_tmp)
     
!     x_prj => SpAMM_tree_2d_symm_times_tree_2d_symm( z%mtx,  x_tmp, z%nxt%tau, &
!                                                     alpha_O=SpAMM_zero, beta_O=SpAMM_one, in_O=x_prj)     
  enddo

end program SpAMM_sandwich_inverse_squareroot
