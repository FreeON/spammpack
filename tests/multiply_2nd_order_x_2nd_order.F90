program test

  use spammpack
  use test_utilities

  integer, parameter :: M = 12
  integer, parameter :: N = 12

  type(spamm_matrix_2nd_order), pointer :: A, B, C
  real(kind(0d0)), dimension(M, N) :: A_dense, C_dense
  integer :: i, j

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => A
  C => spamm_zero_matrix(size(A_dense, 1), size(A_dense, 2))

  call multiply(A, B, C)

  C_dense = matmul(A_dense, A_dense)

  do i = 1, size(C_dense, 1)
    do j = 1, size(C_dense, 2)
      if(abs(C_dense(i, j)-get(C, i, j)) > 1d-10) then
        write(*, *) "matrix element mismatch"
        write(*, "(A,I3,A,I3,A,F7.4)") "C_reference(", i, ",", j, ") = ", C_dense(i, j)
        write(*, "(A,I3,A,I3,A,F7.4)") "          C(", i, ",", j, ") = ", get(C, i, j)
        error stop
      endif
    enddo
  enddo

end program test
