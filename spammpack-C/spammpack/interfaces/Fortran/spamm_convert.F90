module spamm_convert

  use spamm_derived

contains

  function spamm_convert_dense_to_spamm (A)

    type(spamm_matrix) :: spamm_convert_dense_to_spamm
    real(spamm_kind), dimension(:, :), intent(in) :: A

  end function spamm_convert_dense_to_spamm

end module spamm_convert
