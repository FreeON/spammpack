module test_utilities

  use spammpack

  implicit none

contains

  subroutine load_matrix (filename, A)

    character(len = *), intent(in) :: filename
    real(SpAMM_DOUBLE), dimension(:, :), allocatable, intent(inout) :: A

    integer :: N
    integer :: i, j
    real(SpAMM_DOUBLE) :: Aij

    if(allocated(A)) then
      deallocate(A)
    endif

    N = 0
    open(unit = 10, file = filename)
    do while(.true.)
      read(10, *, end = 1) i, j, Aij
      if(i > N) then
        N = i
      endif
      if(j > N) then
        N = j
      endif
    enddo
1   continue

    if(N == 0) then
      write(*, *) "no matrix elements found in matrix file"
    else
      allocate(A(N, N))
      rewind(10)
      do while(.true.)
        read(10, *, end = 2) i, j, Aij
        A(i, j) = Aij
      enddo
2     close(10)
    endif

  end subroutine load_matrix

  subroutine print_matrix (A)

    integer :: i, j
    real(SpAMM_KIND), dimension(:, :), allocatable, intent(in) :: A

    do i = 1, size(A, 1)
      do j = 1, size(A, 2)
        write(*, *) i, j, A(i, j)
      enddo
    enddo

  end subroutine print_matrix

end module test_utilities
