program test

  use spammpack
  use test_utilities

#include <spamm_utility_macros.h>

  implicit none

  integer, parameter :: N = 129

  type(spamm_matrix_order_2), pointer :: A, B, C
  real(kind(0d0)), dimension(N, N) :: A_dense, C_dense
  real(kind(0d0)) :: reference_complexity
  integer :: i, j

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => A
  C => spamm_zero_matrix(size(A_dense, 1), size(A_dense, 2))

  call multiply(A, B, C)
  reference_complexity = dble(N+SPAMM_BLOCK_SIZE-mod(N, SPAMM_BLOCK_SIZE))**3

  write(*, *) "complexity      = "//to_string(C%number_operations)
  write(*, *) "ref. complexity = "//to_string(reference_complexity)
  write(*, *) "complexity/N^3  = "//to_string(C%number_operations/dble(N)**3)

  if(abs(C%number_operations-reference_complexity) > 1e-10) then
    LOG_FATAL("complexity count wrong")
  endif

  C_dense = matmul(A_dense, A_dense)

  do i = 1, size(C_dense, 1)
    do j = 1, size(C_dense, 2)
      if(abs(C_dense(i, j)-get(C, i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("C_reference("//to_string(i)//","//to_string(j)//") = "//to_string(C_dense(i, j)))
        LOG_FATAL("          C("//to_string(i)//","//to_string(j)//") = "//to_string(get(C, i, j)))
        error stop
      endif
    enddo
  enddo

end program test
