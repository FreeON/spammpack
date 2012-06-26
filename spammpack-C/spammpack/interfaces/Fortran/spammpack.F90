module spammpack

  use spamm_derived

  IMPLICIT NONE

  !> Interface to spamm_version().
  interface spamm_version
    function spamm_version ()
    end function spamm_version
  end interface spamm_version

  !> Interface to spamm_exit().
  interface spamm_exit
    subroutine spamm_exit (exitcode)
      integer :: exitcode
    end subroutine spamm_exit
  end interface spamm_exit

  !> Interface to spamm_convert_dense_to_spamm_interface().
  interface spamm_convert_dense_to_spamm_interface
    subroutine spamm_convert_dense_to_spamm_interface (A, B)
      use spamm_derived
      real(SpAMM_KIND), dimension(:, :), intent(in) :: A
      type(SpAMM_Matrix), intent(inout) :: B
    end subroutine spamm_convert_dense_to_spamm_interface
  end interface spamm_convert_dense_to_spamm_interface

  interface add
    module procedure add_spamm_spamm
  end interface add

  interface multiply
    module procedure multiply_spamm_spamm
  end interface multiply

  interface new
    module procedure new_spamm
  end interface new

contains

  type(SpAMM_Matrix) function spamm_convert_dense_to_spamm (A) result (B)
    real(SpAMM_KIND), dimension(:, :), intent(in) :: A
    call spamm_convert_dense_to_spamm_interface(A, B)
  end function spamm_convert_dense_to_spamm

  subroutine multiply_spamm_spamm (A, B, C, tolerance)
    type(SpAMM_Matrix), intent(in) :: A
    type(SpAMM_Matrix), intent(in) :: B
    type(SpAMM_Matrix), intent(inout) :: C
    real(SpAMM_KIND), intent(in), optional :: tolerance
  end subroutine multiply_spamm_spamm

  subroutine add_spamm_spamm (A, B, alpha, beta)
    type(SpAMM_Matrix), intent(in) :: A
    type(SpAMM_Matrix), intent(inout) :: B
    real(SpAMM_KIND), intent(in) :: alpha
    real(SpAMM_KIND), intent(in) :: beta
  end subroutine add_spamm_spamm

  subroutine new_spamm (A)
    type(SpAMM_Matrix), intent(inout) :: A
    call spamm_new_interface(A)
  end subroutine new_spamm

end module spammpack
