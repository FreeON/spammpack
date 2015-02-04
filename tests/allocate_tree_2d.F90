module a_mod
contains

  subroutine test (A, B)

    use spammpack

    type(tree_2d), pointer, intent(in) :: A
    type(tree_2d), pointer, intent(out) :: B

    write(*, "(A)") "entering test()"

    B => new_tree_2d(A%decoration%N)

    write(*, "(A)") "leaving test()"

  end subroutine test

end module a_mod

program allocate_tree_1d

  use a_mod
  use spammpack

  implicit none

  integer, parameter :: N = 10
  type(tree_2d), pointer :: A, B

  A => new_tree_2d(N)
  nullify(B)

  write(*, "(A)") "calling test()"
  call test(A, B)

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

end program allocate_tree_1d
