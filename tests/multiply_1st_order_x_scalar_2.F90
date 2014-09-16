program test

  use spammpack

  implicit none

  integer, parameter :: N = 102
  real(kind(0d0)), parameter :: alpha = 1.2

  type(spamm_matrix_order_1), pointer :: V => null()
  real(kind(0d0)), dimension(N/10) :: V_temp
  real(kind(0d0)), dimension(N) :: V_dense
  integer :: i

  call random_number(V_temp)
  V_dense = 0
  V_dense(1:size(V_temp)) = V_temp

  V => spamm_convert_dense_to_order_1(V_dense)

  call multiply(V, alpha)
  V_dense = alpha*V_dense

  do i = 1, N
    if(abs(V_dense(i)-get(V, i)) > 1d-10) then
      LOG_FATAL("vector element mismatch")
      LOG_FATAL("V_reference("//to_string(i)//") = "//to_string(V_dense(i)))
      LOG_FATAL("          V("//to_string(i)//") = "//to_string(get(V, i)))
      error stop
    endif
  enddo

end program test
