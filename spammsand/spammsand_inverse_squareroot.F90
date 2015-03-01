! Nested SpAMM solvers (SpAMM sandwitches) for matrix functions
program SpAMMSand_nested_inverse_squareroot

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
  allocate(z%tau2(1:slices+1))
  allocate(z%trix(1:slices))
  do i=1,slices+1
     z%tau(i) =  10d0**( logtau_strt + logtau_dlta*float(i-1) )
     z%trix(i) => SpAMM_new_top_tree_1d(s%frill%ndimn)
  enddo

  x_nxt => SpAMM_tree_2d_symm_copy_tree_2d_symm( s, x_in)
  x_tmp => SpAMM_new_top_tree_1d(s%frill%ndimn)

  do i=1,slices

     call spammsand_scaled_newton_shulz_inverse_squareroot( x_nxt, z%trix(i), z%tau(i) )

     x_tmp=>SpAMM_tree_2d_symm_times_tree_2d_symm( x_nxt,      z%trix(i), z%tau(i+1), &
          alpha_O=, beta_O=, c_O=x_tmp)
     
     x_nxt=>SpAMM_tree_2d_symm_times_tree_2d_symm( z%trix(i),  x_tmp,     z%tau(i+1), &
          alpha_O=, beta_O=, c_O=x_nxt)
     
  enddo

end program SpAMMSand_inverse_squareroot
