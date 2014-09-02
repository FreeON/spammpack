program test

  use spammpack

  implicit none

  integer, parameter :: M = 14
  integer, parameter :: N = 14
  real(kind(0d0)), parameter :: alpha = 1.2

  type(spamm_matrix_2nd_order), pointer :: A
  real(kind(0d0)), dimension(M, N) :: A_dense
  integer :: i, j

  call random_number(A_dense)

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)

  do i = 1, M
    A_dense(i, i) = A_dense(i, i)+alpha
  enddo

  call add(A, alpha)

  do i = 1, size(A_dense, 1)
    do j = 1, size(A_dense, 2)
      if(abs(A_dense(i, j)-get(A, i, j)) > 1d-10) then
        call write_log(FATAL, "matrix element mismatch, "// &
          "A_reference("//to_string(i)//","//to_string(j)//") = "//to_string(A_dense(i, j))//", "// &
          "          A("//to_string(i)//","//to_string(j)//") = "//to_string(get(A, i, j)))
      endif
    enddo
  enddo

end program test
