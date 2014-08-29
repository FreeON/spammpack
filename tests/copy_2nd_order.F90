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
        call write_log(FATAL, [ "A("//to_string(i)//","//to_string(j)//") mismatch", &
          "           found "//to_string(Bij), &
          "should have been "//to_string(Aij) ])
      endif
    enddo
  enddo

  call delete(A)
  call delete(B)

end program test
