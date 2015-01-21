module a_mod
contains

  subroutine test (A, B)

    use spammpack

    type(spamm_matrix_1d), intent(in) :: A
    type(spamm_matrix_1d), intent(out) :: B

#ifdef HAVE_CONSTRUCTOR
    B = spamm_matrix_1d(A%N)
#else
    B = new_matrix_1d(A%N)
#endif

  end subroutine test

end module a_mod

program allocate_matrix_1d

  use a_mod
  use spammpack

  implicit none

  type(spamm_matrix_1d) :: A, B

  write(*, "(A)") "calling test()"
  call test(A, B)

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

end program allocate_matrix_1d
