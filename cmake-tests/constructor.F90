module a_mod

  implicit none

  type :: a_t
     integer :: i
  end type a_t

  interface a_t
     module procedure new_a_t
  end interface a_t

contains

  type(a_t) function new_a_t (i) result(a)

    integer, intent(in) :: i

    a%i = i

  end function new_a_t

end module a_mod

program test

  use a_mod

  implicit none

  type(a_t) :: a

  a = a_t(1)

end program test
