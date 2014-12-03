program spamm_trace

  use, intrinsic :: iso_c_binding
  use spammpack

  type(c_ptr) :: A

  real*8, allocatable, dimension(:,:) :: ADense
  integer, dimension(2) :: N = (/ 513, 513 /)
  integer :: chunkTier = 2
  logical :: useLinearTree = .false.
  integer :: i
  real*8 :: trace_reference

  allocate(ADense(N(1), N(2)))
  call random_number(ADense)

  call SpAMM_SetEq(A, ADense, chunkTier, useLinearTree)

  trace_reference = 0D0

  do i = 1, N(1)
    trace_reference = trace_reference+ADense(i,i)
  enddo

end program spamm_trace
