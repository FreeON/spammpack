module test_utilities

  implicit none

  interface mmio

    subroutine mmread(iunit, rep, field, symm, rows, cols, nnz, nnzmax, indx, jndx,ival,rval,cval)
      integer, intent(in) :: iunit
      character(len = 10), intent(inout) :: rep
      character(len = 7), intent(inout) :: field
      character(len = 19), intent(inout) :: symm
      integer, intent(inout) :: rows, cols, nnz
      integer, intent(in) :: nnzmax
      integer, dimension(*), intent(inout) :: indx, jndx, ival
      real(kind(0d0)), dimension(*), intent(inout) :: rval
      complex(kind(0d0)), dimension(*), intent(inout) :: cval
    end subroutine mmread

    subroutine mminfo(iunit,rep,field,symm,rows,cols,nnz)
      integer, intent(in) :: iunit
      character(len = 10), intent(inout) :: rep
      character(len = 7), intent(inout) :: field
      character(len = 19), intent(inout) :: symm
      integer, intent(inout) :: rows, cols, nnz
    end subroutine mminfo

  end interface mmio

contains

  !> Read a matrix in MatrixMarket format.
  !!
  !! @param filename The filename.
  !! @param A The dense matrix
  subroutine read_MM (filename, A)

    implicit none

    character(len = *), intent(in) :: filename
    real(kind(0d0)), allocatable, dimension(:, :), intent(inout) :: A

    integer :: M
    integer :: N
    integer :: number_nonzero, max_number_nonzero
    integer :: i

    character(len = 10) :: rep
    character(len = 7) :: field
    character(len = 19) :: symm
    integer, allocatable, dimension(:) :: indx, jndx
    integer, allocatable, dimension(:) :: ival
    real(kind(0d0)), allocatable, dimension(:) :: rval
    complex(kind(0d0)), allocatable, dimension(:) :: cval

    write(*, "(A)") "reading matrix from "//trim(filename)

    open(unit = 20, file = filename)
    call mminfo(20, rep, field, symm, M, N, number_nonzero)

    write(*, *) "M       = ", M
    write(*, *) "N       = ", N
    write(*, *) "nnzero  = ", number_nonzero
    write(*, *) "%fillin = ", number_nonzero/(dble(M)*dble(N))
    write(*, *) "rep     = "//trim(rep)
    write(*, *) "field   = "//trim(field)
    write(*, *) "symm    = "//trim(symm)
    write(*, *) "mem(A)  = ", 8*M*N/1024.**2, " MiB"

    allocate(indx(number_nonzero))
    allocate(jndx(number_nonzero))
    allocate(ival(number_nonzero))
    allocate(rval(number_nonzero))
    allocate(cval(number_nonzero))

    max_number_nonzero = number_nonzero

    call mmread(20, rep, field, symm, M, N, number_nonzero, max_number_nonzero, indx, jndx, ival, rval, cval)

    write(*, *) "allocating A"
    allocate(A(M, N))

    write(*, *) "converting formats"
    do i = 1, number_nonzero
      A(indx(i), jndx(i)) = rval(i)
    enddo

    deallocate(ival)
    deallocate(rval)
    deallocate(cval)
    deallocate(jndx)
    deallocate(indx)

  end subroutine read_MM

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
      error stop
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
    real(kind(0d0)), dimension(:, :), allocatable, intent(inout) :: A

    integer :: N
    integer :: i, j
    real(kind(0d0)) :: Aij

    if(allocated(A)) then
      deallocate(A)
    endif

    write(*, "(A,A)") "reading matrix from ", trim(filename)

    N = 0
    open(unit = 1, file = trim(filename), form = "unformatted", recl = 8)
    read(1) Aij
    N = int(Aij)

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

  !> Print a matrix.
  !!
  !! @param A The matrix.
  subroutine print_matrix (A)

    integer :: i, j
    character(len = 20) :: format_string
    real(kind(0d0)), dimension(:, :), intent(in) :: A

    write(format_string, "(A,I3,A)") "(", size(A, 2), "ES10.3)"
    do i = 1, size(A, 1)
      write(*, format_string) (A(i, j), j = 1, size(A, 2))
    enddo

  end subroutine print_matrix

  !> Print a matrix, python style.
  !!
  !! @param A The matrix.
  !! @param variable_name The name of the python variable.
  subroutine print_matrix_python (A, variable_name)

    real(kind(0d0)), dimension(:, :), intent(in) :: A
    character(len = *), intent(in) :: variable_name

    integer :: i, j
    character(len = 2000) :: line_string

    write(*, "(A,A)") variable_name, " = numpy.matrix(["
    do i = 1, size(A, 1)
      write(line_string, "(A)") "  ["
      do j = 1, size(A, 2)
        write(line_string, "(A,ES10.3)") trim(line_string), A(i, j)
        if(j < size(A, 2)) then
          write(line_string, "(A,A)") trim(line_string), ", "
        endif
      enddo
      write(line_string, "(A,A)") trim(line_string), "]"
      if(i < size(A, 1)) then
        write(line_string, "(A,A)") trim(line_string), ","
      endif
      write(*, "(A)") trim(line_string)
    enddo
    write(*, "(A)") "  ])"

  end subroutine print_matrix_python

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
