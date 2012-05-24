module load_matrix

  implicit none

contains

  subroutine load (filename, A)

    character(len = *), intent(in)                             :: filename
    real(kind = 8), dimension(:,:), allocatable, intent(inout) :: A

  end subroutine load

end module load_matrix
