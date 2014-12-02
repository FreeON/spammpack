!> A matmul that counts the number of non-zero products.
!!
!! \f$ C \leftarrow A B \f$
!!
!! @param A Matrix A.
!! @param B Matrix B.
!! @param C Matrix C.
!! @param number_operations The operations count.
subroutine counting_matmul (A, B, C, number_operations)

  use spamm_globals

  implicit none

  real(kind(0d0)), intent(in) :: A(:, :), B(:, :)
  real(kind(0d0)), intent(inout) :: C(:, :)
  real(kind(0d0)), intent(out) :: number_operations

  real(kind(0d0)), allocatable :: A_padded(:, :), B_padded(:, :), C_padded(:, :)
  real(kind(0d0)), allocatable :: A_block(:, :), B_block(:, :)
  real(kind(0d0)) :: norm_A, norm_B

  integer :: N, N_padded
  integer :: i, j, k, i_block, j_block

  N = size(A, 1)
  if(mod(N, SPAMM_BLOCK_SIZE) == 0) then
    N_padded = N
  else
    N_padded = N+(SPAMM_BLOCK_SIZE-mod(N, SPAMM_BLOCK_SIZE))
  endif

  write(*, *) "N_padded", N_padded

  allocate(A_padded(N_padded, N_padded))
  allocate(B_padded(N_padded, N_padded))
  allocate(C_padded(N_padded, N_padded))

  allocate(A_block(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
  allocate(B_block(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))

  A_padded = 0d0
  B_padded = 0d0
  C_padded = 0d0

  A_padded(1:N, 1:N) = A
  B_padded(1:N, 1:N) = B

  number_operations = 0

  do i = 1, N/SPAMM_BLOCK_SIZE
    do j = 1, N/SPAMM_BLOCK_SIZE
      do k = 1, N/SPAMM_BLOCK_SIZE
        A_block = A_padded((i-1)*SPAMM_BLOCK_SIZE+1:i*SPAMM_BLOCK_SIZE, &
          (k-1)*SPAMM_BLOCK_SIZE+1:k*SPAMM_BLOCK_SIZE)
        B_block = B_padded((k-1)*SPAMM_BLOCK_SIZE+1:k*SPAMM_BLOCK_SIZE, &
          (j-1)*SPAMM_BLOCK_SIZE+1:j*SPAMM_BLOCK_SIZE)

        norm_A = 0
        norm_B = 0

        do i_block = 1, SPAMM_BLOCK_SIZE
          do j_block = 1, SPAMM_BLOCK_SIZE
            norm_A = norm_A+A_block(i_block, j_block)**2
            norm_B = norm_B+B_block(i_block, j_block)**2
          enddo
        enddo

        if(sqrt(norm_A*norm_B) > 0) then
          C((i-1)*SPAMM_BLOCK_SIZE+1:i*SPAMM_BLOCK_SIZE, &
            (j-1)*SPAMM_BLOCK_SIZE+1:j*SPAMM_BLOCK_SIZE) = &
            C((i-1)*SPAMM_BLOCK_SIZE+1:i*SPAMM_BLOCK_SIZE, &
            (j-1)*SPAMM_BLOCK_SIZE+1:j*SPAMM_BLOCK_SIZE) + &
            matmul(A_block, B_block)
          number_operations = number_operations+SPAMM_BLOCK_SIZE**3
        endif
      enddo
    enddo
  enddo

  deallocate(A_block)
  deallocate(B_block)

  deallocate(A_padded)
  deallocate(B_padded)
  deallocate(C_padded)

end subroutine counting_matmul

program test

  use spammpack
  use test_utilities

#include <spamm_utility_macros.h>

  implicit none

  interface
    subroutine counting_matmul (A, B, C, number_operations)
      real(kind(0d0)), intent(in) :: A(:, :), B(:, :)
      real(kind(0d0)), intent(inout) :: C(:, :)
      real(kind(0d0)), intent(out) :: number_operations
    end subroutine counting_matmul
  end interface

  logical, parameter :: PRINT_MATRICES = .false.
  integer, parameter :: N = 128
  integer, parameter :: BANDWIDTH = 4
  integer, parameter :: ITERATIONS = 5

#ifdef SPAMM_COUNTERS
  type(spamm_matrix_2nd_order), pointer :: A, B, C
  real(kind(0d0)), dimension(N, N) :: A_dense, C_dense
  real(kind(0d0)) :: reference_norm
  real(kind(0d0)), dimension(ITERATIONS) :: reference_complexity
  integer :: i, j
#endif

#ifndef SPAMM_COUNTERS 
  write(*, "(A)") "skipping test, library was not compiled with SPAMM_COUNTERS"
#else
  A_dense = 0
  do i = 1, N
    do j = 1, N
      if(abs(i-j) <= BANDWIDTH) then
        A_dense(i, j) = 1/dble(N)**(1/3.)
      endif
    enddo
  enddo

  if(PRINT_MATRICES) then
    call print_matrix(A_dense)
  endif

  A => spamm_convert_dense_to_matrix_2nd_order(A_dense)
  B => A
  C => spamm_zero_matrix(N, N)

  write(*, "(A)") "number non-zero = "//to_string(A%number_nonzeros)
  write(*, "(A)") "ref. non-zeros  = "//to_string(count_nonzero(A_dense))

  if(abs(A%number_nonzeros-count_nonzero(A_dense)) > 1d-10) then
    LOG_FATAL("non-zero count wrong")
    error stop
  endif

  do i = 1, ITERATIONS
    ! Compute the reference.
    call counting_matmul(A_dense, A_dense, C_dense, reference_complexity(i))
    !C_dense = matmul(A_dense, A_dense)
    A_dense = C_dense

    ! Compute with SpAMM.
    call multiply(A, B, C)

    write(*, "(A)") to_string(i)//": non-zeros       = "//to_string(C%number_nonzeros)
    write(*, "(A)") to_string(i)//": ref. non-zeros  = "//to_string(count_nonzero(C_dense))
    write(*, "(A)") to_string(i)//": complexity      = "//to_string(C%number_operations)
    write(*, "(A)") to_string(i)//": ref. complexity = "//to_string(reference_complexity(i))
    write(*, "(A)") to_string(i)//": complexity/N^3  = "//to_string(C%number_operations/dble(N)**3)

    if(abs(C%number_nonzeros-count_nonzero(C_dense)) > 1e-10) then
      LOG_FATAL(to_string(i)//": non-zero count wrong")
      error stop
    endif

    if(abs(C%number_operations-reference_complexity(i)) > 1e-10) then
      LOG_FATAL(to_string(i)//": complexity count wrong")
      error stop
    endif

    call copy(C, A)
  enddo

  if(PRINT_MATRICES) then
    call print_matrix(C_dense)
  endif

  reference_norm = matrix_norm(C_dense)
  write(*, "(A)") "norm_ref = "//to_string(reference_norm)
  write(*, "(A)") "norm     = "//to_string(C%norm)
  if(abs(reference_norm-C%norm) > 1d-10) then
    write(*, "(A)") "norm mismatch"
  endif

  do i = 1, size(C_dense, 1)
    do j = 1, size(C_dense, 2)
      if(abs((C_dense(i, j)-get(C, i, j))/C_dense(i, j)) > 1d-10) then
        LOG_FATAL("matrix element mismatch")
        LOG_FATAL("C_reference("//to_string(i)//","//to_string(j)//") = "//to_string(C_dense(i, j)))
        LOG_FATAL("          C("//to_string(i)//","//to_string(j)//") = "//to_string(get(C, i, j)))
        error stop
      endif
    enddo
  enddo

  write(*, *) "done"
#endif

end program test
