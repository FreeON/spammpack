program SetEq

  use, intrinsic :: iso_c_binding
  use spammpack

  type(c_ptr) :: A
  real*8, allocatable, dimension(:,:) :: ADense
  integer, dimension(2) :: N = (/ 513, 513 /)
  integer :: chunkTier = 2
  logical :: useLinearTree = .true.
  integer :: i, j

  allocate(ADense(N(1), N(2)))
  call random_number(ADense)

  call SpAMM_SetEq(A, ADense, chunkTier, useLinearTree)

  do i = 1, N(1)
    do j = 1, N(2)
      if(abs(SpAMM_Get(i, j, A)-ADense(i, j)) > 1e-6) then
        write(*,*) "mismatch"
        write(*,*) "i = ", i
        write(*,*) "j = ", j
        write(*,*) "ADense = ", ADense(i, j)
        write(*,*) "A      = ", SpAMM_Get(i, j, A)
        stop
      endif
    enddo
  enddo

end program SetEq
