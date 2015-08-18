program test

  use spammpack

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: M = 14
  integer, parameter :: N = 14
  real(kind(0d0)), parameter :: alpha = 1.2, beta = 0.3

  type(spamm_matrix_order_2), pointer :: A
  type(spamm_matrix_order_2), pointer :: B
  real(kind(0d0)), dimension(M, N) :: A_dense
  real(kind(0d0)), dimension(M, N) :: B_dense
  integer :: i, j

  call random_number(A_dense)
  call random_number(B_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => spamm_convert_dense_to_matrix_2nd_order(B_dense)

  forall(i = 1:N, j = 1:N)
    A_dense(i, j) = alpha*A_dense(i, j)+beta*B_dense(i, j)
  end forall

  call add(A, B, alpha, beta)

  do i = 1, M
    do j = 1, N
      if(abs(A_dense(i, j)-get(A, i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("A_reference("//to_string(i)//","//to_string(j)//") = "//to_string(A_dense(i, j)))
        LOG_FATAL("          A("//to_string(i)//","//to_string(j)//") = "//to_string(get(A, i, j)))
        error stop
      endif
    enddo
  enddo

  LOG_INFO("passed dense+dense test")

  call delete(A)
  call delete(B)

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => spamm_identity_matrix(M, N)

  do i = 1, M
    do j = 1, N
      A_dense(i, j) = alpha*A_dense(i, j)
      if(i == j) then
        A_dense(i, j) = A_dense(i, j)+beta
      endif
    enddo
  enddo

  call add(A, B, alpha, beta)

  do i = 1, M
    do j = 1, N
      if(abs(A_dense(i, j)-get(A, i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("A_reference("//to_string(i)//","//to_string(j)//") = "//to_string(A_dense(i, j)))
        LOG_FATAL("          A("//to_string(i)//","//to_string(j)//") = "//to_string(get(A, i, j)))
        error stop
      endif
    enddo
  enddo

  LOG_INFO("passed dense+I test")

end program test
