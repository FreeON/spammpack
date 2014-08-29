program test

  use spammpack
  use test_utilities

  implicit none

  integer, parameter :: M = 6
  integer, parameter :: N = 11

  type(spamm_matrix_2nd_order), pointer :: A
  real(kind(0d0)), dimension(M, N) :: A_dense
  real(kind(0d0)) :: Aij
  integer :: i, j

  call random_number(A_dense)

  call print_matrix(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  do i = 1, M
    do j = 1, N
      Aij = get(A, i, j)
      if(abs(Aij-A_dense(i, j)) /= 0) then
        call write_log(FATAL, [ "A("//to_string(i)//","//to_string(j)//") mismatch", &
          "           found "//to_string(Aij), &
          "should have been "//to_string(A_dense(i, j)) ])
      endif
    enddo
  enddo

  call delete(A)

end program test
