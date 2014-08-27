module test_utilities

  implicit none

contains

  subroutine load_matrix (filename, A)

    character(len = *), intent(in) :: filename
    real(kind(0d0)), dimension(:, :), allocatable, intent(inout) :: A

    integer :: N
    integer :: i, j
    real(kind(0d0)) :: Aij

    if(allocated(A)) then
      deallocate(A)
    endif

    write(*, "(A,A)") "reading matrix from ", trim(filename)

    N = 0
    open(unit = 10, file = trim(filename))
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
      write(*, "(A)") "no matrix elements found in matrix file"
      stop
    else
      allocate(A(N, N))
      A = 0.0D0
      rewind(10)
      do while(.true.)
        read(10, *, end = 2) i, j, Aij
        A(i, j) = Aij
      enddo
2     close(10)
    endif

  end subroutine load_matrix

  subroutine load_matrix_binary (filename, A)

    character(len = *), intent(in) :: filename
    real*8, dimension(:, :), allocatable, intent(inout) :: A

    integer :: N
    integer :: i, j
    real*8 :: Aij

    if(allocated(A)) then
      deallocate(A)
    endif

    write(*, "(A,A)") "reading matrix from ", trim(filename)

    N = 0
    open(unit = 1, file = trim(filename), form = "unformatted", recl = 8)
    read(1) Aij
    N = Aij

    !write(*, *) "N = ", N

    allocate(A(N, N))
    do i = 1, N
      do j = 1, N
        read(1) A(i, j)
        !write(*, *) "Aij = ", Aij
      enddo
    enddo

    close(1)

  end subroutine load_matrix_binary

  subroutine print_matrix (A)

    integer :: i, j
    character(len = 20) :: format_string
    real(kind(0d0)), dimension(:, :), intent(in) :: A

    write(format_string, "(A,I3,A)") "(", size(A, 2), "ES10.3)"
    do i = 1, size(A, 1)
      write(*, format_string) (A(i, j), j = 1, size(A, 2))
    enddo

  end subroutine print_matrix

  !> Calculate the Frobenius norm of a dense matrix.
  !!
  !! @param A The matrix
  !!
  !! @return The Frobenius norm.
  function matrix_norm (A) result(norm)

    real(kind(0d0)), dimension(:, :) :: A
    real(kind(0d0)) :: norm
    integer :: i, j

    if(size(A, 1) == size(A, 2)) then
      norm = sqrt(sum(matmul(A, transpose(A)), &
        reshape((/ ((i == j, i = 1, size(A, 1)), j = 1, size(A, 2)) /), &
        (/ size(A, 1), size(A, 2) /))))
    else
      write(*, *) "This implementation can only handle square matrices"
      norm = -1
    endif

  end function matrix_norm

end module test_utilities
