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

  interface add
    module procedure add_spamm_spamm
  end interface add

  interface multiply
    module procedure multiply_spamm_spamm
  end interface multiply

  interface new
    module procedure new_spamm
  end interface new

  interface norm
    module procedure spamm_norm
  end interface norm

  interface spamm_get_norm_interface
    function spamm_get_norm_interface (A) result(norm_F)
      use spamm_derived
      real(SpAMM_SINGLE) :: norm_F
      type(SpAMM_Matrix), intent(in) :: A
    end function spamm_get_norm_interface
  end interface spamm_get_norm_interface

contains

  type(SpAMM_Matrix) function spamm_convert_dense_to_spamm (A) result (B)
    real(SpAMM_SINGLE), dimension(:, :), intent(in) :: A
    call spamm_convert_dense_to_spamm_interface(size(A, 1), size(A, 2), A(1, 1), B)
  end function spamm_convert_dense_to_spamm

  subroutine multiply_spamm_spamm (A, B, C, tolerance)
    type(SpAMM_Matrix), intent(in) :: A
    type(SpAMM_Matrix), intent(in) :: B
    type(SpAMM_Matrix), intent(inout) :: C
    real(SpAMM_SINGLE), intent(in), optional :: tolerance
    call spamm_multiply_spamm_spamm_interface(A, B, C, tolerance)
  end subroutine multiply_spamm_spamm

  subroutine add_spamm_spamm (A, B, alpha, beta)
    type(SpAMM_Matrix), intent(in) :: A
    type(SpAMM_Matrix), intent(inout) :: B
    real(SpAMM_SINGLE), intent(in) :: alpha
    real(SpAMM_SINGLE), intent(in) :: beta
    call spamm_add_spamm_spamm_interface(A, B, alpha, beta)
  end subroutine add_spamm_spamm

  subroutine new_spamm (M, N, A)
    integer, intent(in) :: M, N
    type(SpAMM_Matrix), intent(inout) :: A
    call spamm_new_interface(M, N, A)
  end subroutine new_spamm

  real(SpAMM_SINGLE) function spamm_norm (A) result(norm_F)
    type(SpAMM_Matrix), intent(in) :: A
    norm_F = spamm_get_norm_interface(A)
  end function spamm_norm

end module spammpack
