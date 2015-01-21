module a_mod
contains

  subroutine test (A, B)

    use spammpack

    type(spamm_tree_1d), intent(in) :: A
    type(spamm_tree_1d), intent(out) :: B

    write(*, "(A)") "entering test()"

#ifdef HAVE_CONSTRUCTOR
    B = spamm_tree_1d(A%decoration%N)
#else
    B = new_tree_1d(A%decoration%N)
#endif

    write(*, "(A)") "leaving test()"

  end subroutine test

end module a_mod

program allocate_tree_1d

  use a_mod
  use spammpack

  implicit none

  integer, parameter :: N = 10
  type(spamm_tree_1d) :: A, B

#ifdef HAVE_CONSTRUCTOR
  A = spamm_tree_1d(N)
#else
  A = new_tree_1d(N)
#endif

  write(*, "(A)") "calling test()"
  call test(A, B)

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

end program allocate_tree_1d
