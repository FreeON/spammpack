module a_mod
contains

  subroutine test (A, B)

    use spammpack

    type(tree_2d_symmetric), pointer, intent(in) :: A
    type(tree_2d_symmetric), pointer, intent(out) :: B

    write(*, "(A)") "entering test()"

    B => new_tree_2d_symmetric(A%decoration%N)

    write(*, "(A)") "leaving test()"

  end subroutine test

end module a_mod

program testprogram

  use a_mod
  use spammpack

  implicit none

  integer, parameter :: N = 10
  type(tree_2d_symmetric), pointer :: A, B

  A => new_tree_2d_symmetric(N)
  nullify(B)

  write(*, "(A)") "calling test()"
  call test(A, B)

  write(*, "(A)") "A: "//trim(A%to_string())
  write(*, "(A)") "B: "//trim(B%to_string())

  call delete_tree_2d_symmetric(A)
  call delete_tree_2d_symmetric(B)

end program testprogram
