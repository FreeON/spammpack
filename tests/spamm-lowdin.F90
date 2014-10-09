program spamm_lowdin

  use spammpack
  use test_utilities

  implicit none

  real(kind(0d0)), parameter :: tolerance = 1d-8

  type(spamm_matrix_2nd_order), pointer :: S => null(), X => null()
  real(kind(0d0)), dimension(:, :), allocatable :: S_dense
  character(len = 1000) :: matrix_filename
  real :: start_time, end_time

  if(command_argument_count() == 0) then
    write(*, "(A)") "Please specify the matrix file to use (the overlap matrix)"
    error stop
  else if(command_argument_count() > 1) then
    write(*, "(A)") "Too many command line arguments. One is plenty."
    error stop
  else
    call get_command_argument(1, matrix_filename)
  endif

  call read_MM(matrix_filename, S_dense)
  S => spamm_convert_dense_to_matrix_2nd_order(S_dense)

  call cpu_time(start_time)
  X => spamm_inverse_sqrt_schulz(S, tolerance)
  call cpu_time(end_time)

  write(*, "(A)") "X (S^{-1/2}) fillin: "//to_string(X%number_nonzeros)
  write(*, "(A)") "CPU time: "//to_string(end_time-start_time)

end program spamm_lowdin
