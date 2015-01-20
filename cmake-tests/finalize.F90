module test_mod

  implicit none

  type :: a_t
   contains
     final :: a_destructor
  end type a_t

contains

  elemental subroutine a_destructor (a)
    type(a_t), intent(inout) :: a
  end subroutine a_destructor

end module test_mod

program test

  use test_mod

  implicit none

  type(a_t) :: a

end program test
