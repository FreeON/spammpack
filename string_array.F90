module testmod

  implicit none

  contains

    subroutine foo (msg)

      character(len = *), dimension(:) :: msg
      integer :: i

      do i = 1, size(msg)
        write(*, *) len(msg(i)), "'"//msg(i)//"'"
      enddo

    end subroutine foo

end module testmod

program test

  use testmod

  call foo([ "0", "01", "012" ])

end program test
