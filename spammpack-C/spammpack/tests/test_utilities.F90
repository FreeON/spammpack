module test_utilities

  use spammpack

  implicit none

  interface load_matrix
    module procedure load_matrix_single, load_matrix_double
  end interface load_matrix

contains

  subroutine load_matrix_single (filename, A)

    character(len = *), intent(in) :: filename
    real(SpAMM_SINGLE), dimension(:, :), allocatable, intent(inout) :: A
    real(SpAMM_DOUBLE), dimension(:, :), allocatable :: A_double
    integer :: i, j

    call load_matrix_double(filename, A_double)

    allocate(A(size(A_double, 1), size(A_double, 2)))
    do i = 1, size(A_double, 1)
      do j = 1, size(A_double, 1)
        A(i, j) = A_double(i, j)
      enddo
    enddo

  end subroutine load_matrix_single

  subroutine load_matrix_double (filename, A)

    character(len = *), intent(in) :: filename
    real(SpAMM_DOUBLE), dimension(:, :), allocatable, intent(inout) :: A

    integer :: N
    integer :: i, j
    real(SpAMM_DOUBLE) :: Aij

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

  end subroutine load_matrix_double

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
    !open(unit = 1, file = trim(filename), form = "unformatted", access = "direct", recl = 8)
    !open(unit = 1, file = trim(filename), form = "binary", recl = 8)
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
    real(SpAMM_KIND), dimension(:, :), allocatable, intent(in) :: A

    do i = 1, size(A, 1)
      do j = 1, size(A, 2)
        write(*, *) i, j, A(i, j)
      enddo
    enddo

  end subroutine print_matrix

end module test_utilities
