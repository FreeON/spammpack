module spamm_convert

  use spamm_derived

contains

  type(SpAMM_Matrix) function spamm_convert_dense_to_spamm (A) result (B)

    real(SpAMM_KIND), dimension(:, :), intent(in) :: A

  end function spamm_convert_dense_to_spamm

end module spamm_convert
