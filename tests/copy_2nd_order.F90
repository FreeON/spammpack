program test

  use spammpack

  implicit none

  integer, parameter :: M = 5
  integer, parameter :: N = 6

  type(spamm_matrix_2nd_order), pointer :: A, B
  real(kind(0d0)), dimension(M, N) :: A_dense
  integer :: i, j
  real(kind(0d0)) :: Aij, Bij

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => null()

  call copy(A, B)

  do i = 1, M
    do j = 1, N
      Aij = get(A, i, j)
      Bij = get(B, i, j)
      if(abs(Aij-Bij) /= 0) then
        write(*, *) "A(", i, ",", j, ") mismatch"
        write(*, *) "           found ", Bij
        write(*, *) "should have been ", Aij
        error stop
      endif
    enddo
  enddo

end program test
