program spamm_lowdin

  use spammpack
  use test_utilities

  implicit none

  type(spamm_matrix_2nd_order), pointer :: S => null()
  real(kind(0d0)), dimension(:, :), allocatable :: S_dense
  character(len = 1000) :: matrix_filename

  if(command_argument_count() == 0) then
    write(*, "(A)") "Please specify the matrix file to use (the overlap matrix)"
  else if(command_argument_count() > 1) then
    write(*, "(A)") "Too many command line arguments. One is plenty."
  else
    call get_command_argument(1, matrix_filename)
  endif

  call read_MM(matrix_filename, S_dense)

end program spamm_lowdin
