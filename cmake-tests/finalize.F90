module test_mod

  implicit none

  type :: a_t
   contains
     final :: a_destructor
  end type a_t

contains

  function a_destructor (a)
    type(a_t) :: a
  end function a_destructor

end module test_mod
