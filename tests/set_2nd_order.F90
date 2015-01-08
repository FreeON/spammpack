program test

  use spammpack

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: N = 5

  type(spamm_matrix_order_2), pointer :: A
  real(kind(0d0)), dimension(N, N) :: A_dense
  integer :: i, j

  call random_number(A_dense)

  call new(N, N, A)

  do i = 1, N
    do j = 1, N
      call set(A, i, j, A_dense(i, j))
    enddo
  enddo

  do i = 1, size(A_dense, 1)
    do j = 1, size(A_dense, 2)
      if(abs(A_dense(i, j)-get(A, i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("A_reference("//to_string(i)//","//to_string(j)//") = "//to_string(A_dense(i, j)))
        LOG_FATAL("          C("//to_string(i)//","//to_string(j)//") = "//to_string(get(A, i, j)))
        error stop
      endif
    enddo
  enddo

end program test
