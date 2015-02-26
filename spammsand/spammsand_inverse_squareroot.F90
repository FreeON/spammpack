! Nested SpAMM solvers (SpAMM sandwitches) for matrix functions

program SpAMMSand_inverse_squareroot

  USE spammpack 
  USE sandpack 
  use test_utilities
  
  implicit none

  TYPE(spammsand_tree_2d_symm) :: z

  real(kind(0d0)), parameter ::   Tau_1=1d-4
  real(kind(0d0)), parameter ::   Tau_2=1d-6
  real(kind(0d0)), parameter ::   Tau_3=1d-8
  real(kind(0d0)), parameter ::   Tau_4=1d-10
  real(kind(0d0)), parameter ::   Tau_5=1d-12
  real(kind(0d0)), parameter ::   Tau_6=1d-14

  type(SpAMM_tree_2d_symm), pointer :: S => null()
  real(SpAMM_KIND),     allocatable :: S_dense(:, :)

  character(len = 1000)             :: matrix_filename
  real(SpAMM_KIND)                  :: x_hi

!  real :: start_time, end_time

  call get_command_argument(1, matrix_filename)

  call read_MM(matrix_filename, S_dense)

  s => SpAMM_convert_dense_to_tree_2d_symm(S_DENSE) 

  CALL SpAMM_print_tree_2d_symm_recur (s) 

  x_hi=SpAMMSand_rqi_extremal(s,1d-4,high_O=.TRUE.)

  WRITE(*,*)' X + = ',x_hi




!     CALL GSOLVE_NORMALIZE_MATRIX(S, EvMax )
!  ENDIF

!!$  CALL Copy(S,X1)
!!$  CALL GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(X1, Z1, Tau_1 )
!!$
!!$  M=S%M
!!$  N=S%N
!!$  ZT => SpAMM_identity_matrix(M,N)
!!$  T1 => SpAMM_identity_matrix(M,N)
!!$  CALL Multiply( ZT, SpAMM_Zero)       
!!$  CALL Multiply( T1, SpAMM_Zero)       
!!$
!!$  CALL SpAMM_Transpose_QuTree( Z1%root, ZT%root )
!!$
!!$  CALL Copy(S,X2)
!!$
!!$  CALL spamm_convert_order_2_to_dense (Z1, Z_dense)       
!!$  ZT=>SpAMM_convert_dense_to_matrix_2nd_order(TRANSPOSE(Z_dense))
!!$
!!$  CALL Multiply( X2,  Z1, T1, Tau_2 )
!!$  CALL Multiply( ZT, T1,  X2, Tau_2      )
!!$
!!$  CALL GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(X2, Z2, Tau_2 )
!!$
!!$  STOP
!!$
!!$111 CONTINUE
!!$
!!$
!!$  evec=S_dense
!!$  CALL SetScalarSpectrum(1d10,eval,N)
!!$  DO i=1,N
!!$     S_dense(:,I)=S_dense(:,i)*eval(i)
!!$  ENDDO
!!$  SNew=MATMUL(evec,TRANSPOSE(S_dense))
!!$  S => SpAMM_convert_dense_to_matrix_2nd_order(SNew)
!!$
!!$  call cpu_time(start_time)
!!$  call GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT(S, Z, 1d-3)!tolerance) 
!!$  call cpu_time(end_time)

  STOP  

end program SpAMMSand_inverse_squareroot
