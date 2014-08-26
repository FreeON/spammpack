program test

  use spammpack
  use test_utilities

  integer, parameter :: N = 12

  integer, parameter :: iterations = 3

  type(spamm_matrix_2nd_order), pointer :: A, B, C
  real(kind(0d0)), dimension(N, N) :: A_dense, C_dense
  real(kind(0d0)), dimension(N) :: diagonal
  integer :: i, j

  call random_number(diagonal)
  A_dense = 0
  do i = 1, N
    A_dense(i, i) = diagonal(i)
  enddo

  ! Tri-diagonal matrix.
  do i = 1, N
    if(i < N) then
      A_dense(i, i+1) = A_dense(i, i)/2.0
    endif
    if(i > 1) then
      A_dense(i, i-1) = A_dense(i, i)/2.0
    endif
  enddo

  call print_matrix(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => A
  C => spamm_zero_matrix(size(A_dense, 1), size(A_dense, 2))

  do i = 1, iterations
    call multiply(A, B, C)
    call copy(C, A)
    C_dense = matmul(A_dense, A_dense)
    A_dense = C_dense
  enddo

  do i = 1, size(C_dense, 1)
    do j = 1, size(C_dense, 2)
      if(abs((C_dense(i, j)-get(C, i, j))/C_dense(i, j)) > 1d-10) then
        write(*, *) "matrix element mismatch"
        write(*, "(A,I3,A,I3,A,ES11.4)") "C_reference(", i, ",", j, ") = ", C_dense(i, j)
        write(*, "(A,I3,A,I3,A,ES11.4)") "          C(", i, ",", j, ") = ", get(C, i, j)
        error stop
      endif
    enddo
  enddo

end program test
