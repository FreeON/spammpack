program matrix_multiply

  use lal

  implicit none

  integer, parameter :: M = 4
  integer, parameter :: N = 4

  integer          :: i, j, k
  type(lal_matrix) :: A, B, C, C_reference

  if(f90_lal_allocate(M, N, A) /= 0) THEN
    call exit(-1)
  endif
  if(f90_lal_allocate(M, N, B) /= 0) THEN
    call exit(-1)
  endif
  if(f90_lal_allocate(M, N, C) /= 0) THEN
    call exit(-1)
  endif
  if(f90_lal_allocate(M, N, C_reference) /= 0) THEN
    call exit(-1)
  endif

  call f90_lal_rand(A)
  call f90_lal_rand(B)

  ! Multiply by hand.
  call f90_lal_zero(C_reference)
  do i = 1, M
    do j = 1, N
      do k = 1, M
        call f90_lal_set(i, j, f90_lal_get(i, j, C_reference) &
          + f90_lal_get(i, k, A)*f90_lal_get(k, j, B), C_reference)
      enddo
    enddo
  enddo

  ! Multiply by library.
  call f90_lal_dgemm("N", "N", M, N, N, 1.0d0, A, N, B, M, 1.0d0, C, N)

  if(f90_lal_equals(C, C_reference) /= 0) then
    write(*,"(A)") "[f90_matrix_multiply] C is not equal to C_reference"
    write(*,"(A)") "[f90_matrix_multiply] C ="
    call f90_lal_print(C)

    write(*,"(A)") "[f90_matrix_multiply] C_reference ="
    call f90_lal_print(C_reference)
    call exit(-1)
  endif

end program matrix_multiply
