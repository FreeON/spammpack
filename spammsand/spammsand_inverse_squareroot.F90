module SpAMMsand_inverse_squareroot

  USE spammpack 

  implicit none

contains

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot( x, z, tau , t)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(INOUT) :: x, z, t
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau

    INTEGER                                           :: i,j, n
    REAL(SpAMM_KIND)                                  :: EvMin, EvMax, sc, xo, xo_prev 
    REAL(SpAMM_KIND)                                  :: xo_analytic, xmax,delta, FillN, FillN_prev

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
  use test_utilities
  
  implicit none

  TYPE(spammsand_tree_2d_symm)      :: z

  type(SpAMM_tree_2d_symm), pointer :: S => null()

  character(len = 1000)             :: matrix_filename

  real(SpAMM_KIND)                  :: x_hi

  integer :: i,j

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)

  ! matrix to inverse factor
  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 

  ! the max eigenvalue
  x_hi=SpAMMSand_rqi_extremal(s,1d-4,high_O=.TRUE.)

  ! normalize the max ev of s to 1.  
  CALL SpAMM_scalar_times_tree_2d_symm(SpAMM_one/x_hi, s)

  logtau_strt=-4
  logtau_stop=-12
  logtau_dlta=(logtau_stop-logtau_strt)/dble(slices+1) !????????????
  
  ! initialize the slices
  allocate(z%tau(1:slices+1))
  allocate(z%trix(1:slices))
  do i=1,slices+1
     z%tau(i) =  10d0**( logtau_strt + logtau_dlta*float(i-1) )
     z%trix(i) => SpAMM_new_top_tree_1d(s%frill%ndimn)
  enddo

  x_nxt => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, x_nxt )
  x_tmp => SpAMM_new_top_tree_1d( s%frill%ndimn)

  do i=1,slices

     call spammsand_scaled_newton_shulz_inverse_squareroot( x_nxt, z%trix(i), z%tau(i), x_tmp )

     x_tmp => SpAMM_tree_2d_symm_times_tree_2d_symm( x_nxt,      z%trix(i), z%tau(i+1), &
                                                     alpha_O=SpAMM_zero, beta_O=SpAMM_one, c_O=x_tmp)
     
     x_nxt => SpAMM_tree_2d_symm_times_tree_2d_symm( z%trix(i),  x_tmp,     z%tau(i+1), &
                                                     alpha_O=SpAMM_zero, beta_O=SpAMM_one, c_O=x_nxt)     
  enddo

end program SpAMM_sandwich_inverse_squareroot
