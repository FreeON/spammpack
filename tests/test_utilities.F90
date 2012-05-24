module test_utilities

  implicit none

contains

  subroutine load_matrix (filename, A)

    character(len = *), intent(in)                             :: filename
    real(kind = 8), dimension(:,:), allocatable, intent(inout) :: A

  end subroutine load_matrix

  subroutine dgemm (N, alpha, A, B, beta, C)

    integer :: i, j, k
    integer, intent(in) :: N
    real(kind = 8), intent(in) :: alpha
    real(kind = 8), dimension(:,:), allocatable, intent(in) :: A
    real(kind = 8), dimension(:,:), allocatable, intent(in) :: B
    real(kind = 8), intent(in) :: beta
    real(kind = 8), dimension(:,:), allocatable, intent(inout) :: C

    do i = 1, N
      do j = 1, N
        C(i, j) = beta*C(i, j)
        do k = 1, N
          C(i, j) = C(i, j) + alpha*A(i, k)*B(k, j)
        enddo
      enddo
    enddo

  end subroutine dgemm

end module test_utilities
