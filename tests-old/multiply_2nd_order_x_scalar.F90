program test

  use spammpack
  use test_utilities

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: M = 12
  integer, parameter :: N = 12
  real(kind(0d0)), parameter :: alpha = 1.2

  type(spamm_matrix_order_2), pointer :: A
  real(kind(0d0)), dimension(M, N) :: A_dense
  integer :: i, j

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  call multiply(A, alpha)
  A_dense = alpha*A_dense

  do i = 1, size(A_dense, 1)
    do j = 1, size(A_dense, 2)
      if(abs(A_dense(i, j)-get(A, i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("A_reference("//to_string(i)//","//to_string(j)//") = "//to_string(A_dense(i, j)))
        LOG_FATAL("          A("//to_string(i)//","//to_string(j)//") = "//to_string(get(A, i, j)))
        error stop
      endif
    enddo
  enddo

end program test
