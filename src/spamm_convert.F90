!> Defines conversion operation between different data structures and SpAMM.
!!
!! @copyright
!!
!! Copyright (c) 2014, Los Alamos National Laboratory
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice, this
!! list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!! this list of conditions and the following disclaimer in the documentation
!! and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its contributors
!! may be used to endorse or promote products derived from this software without
!! specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! @author Matt Challacombe matt.challacombe@freeon.org
!! @author Nicolas Bock nicolas.bock@freeon.org
MODULE SpAMM_CONVERT

  use spamm_globals
  use spamm_types
  USE SpAMM_ALGEBRA
  use spamm_management

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: spamm_convert_dense_to_matrix_2nd_order

  CONTAINS

  !> Recursively convert a dense matrix to a quadtree.
  !!
  !! @param A The dense matrix.
  !! @param qA A pointer to a quadtree node.
  !! @param i_lower The lower value of the row index.
  !! @param i_upper The upper value of the row index.
  !! @param j_lower The lower value of the column index.
  !! @param j_upper The upper value of the column index.
  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree (A, qA, &
      i_lower, i_upper, j_lower, j_upper)

    REAL(SpAMM_KIND), DIMENSION(:,:), INTENT(IN) :: A
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    INTEGER, INTENT(IN) :: i_lower, i_upper
    INTEGER, INTENT(IN) :: j_lower, j_upper

    INTEGER :: i, j, i_dense, j_dense
    INTEGER :: A_rows, A_cols
    integer :: convert_rows, convert_columns

    A_rows = i_upper-i_lower+1
    A_cols = j_upper-j_lower+1

    if(.not. associated(qA)) then
      allocate(qA)
      qA%i_lower = i_lower
      qA%i_upper = i_upper
      qA%j_lower = j_lower
      qA%j_upper = j_upper
    endif

    !write(*, *) "q: ", i_lower, i_upper, j_lower, j_upper

    IF(A_rows <= SPAMM_BLOCK_SIZE .AND. A_cols <= SPAMM_BLOCK_SIZE)THEN
      IF(A_rows < SPAMM_BLOCK_SIZE .OR. A_cols < SPAMM_BLOCK_SIZE) THEN
        WRITE(*, *) "LOGIC ERROR IN SpAMM: padding error"
        WRITE(*, *) "A_rows = ", A_rows
        WRITE(*, *) "A_cols = ", A_cols
        WRITE(*, *) "SPAMM_BLOCK_SIZE = ", SPAMM_BLOCK_SIZE
        error stop
      ELSE
        if(allocated(qA%blok)) then
          deallocate(qA%blok)
        endif

        ALLOCATE(qA%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))

        ! Set new block to zero.
        qA%Blok = 0

        ! Copy matrix elements.
        convert_rows = min(i_upper, size(A, 1))
        convert_columns = min(j_upper, size(A, 2))

        qA%Blok(1:convert_rows-i_lower+1, 1:convert_columns-j_lower+1) = &
          A(i_lower:convert_rows, j_lower:convert_columns)

        !write(*, *) "q%blok"
        !do i = 1, SPAMM_BLOCK_SIZE
        !  write(*, "(4E10.3)") (qA%blok(i, j), j = 1, SPAMM_BLOCK_SIZE)
        !enddo

        ! Count number non-zeros.
        DO i = 1, SPAMM_BLOCK_SIZE
          DO j = 1, SPAMM_BLOCK_SIZE
            IF(qA%Blok(i, j) /= 0.0D0) THEN
              qA%number_nonzeros = qA%number_nonzeros+1
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ELSE
      ! Avoid slicing here for performance.
      CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad11, &
        i_lower, &
        i_lower+A_rows/2-1, &
        j_lower, &
        j_lower+A_cols/2-1)
      CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad12, &
        i_lower, &
        i_lower+A_rows/2-1, &
        j_lower+A_cols/2, &
        j_upper)
      CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad21, &
        i_lower+A_rows/2, &
        i_upper, &
        j_lower, &
        j_lower+A_cols/2-1)
      CALL SpAMM_Convert_Dense_2_QuTree(A, qA%Quad22, &
        i_lower+A_rows/2, &
        i_upper, &
        j_lower+A_cols/2, &
        j_upper)

      qA%number_nonzeros = &
        qA%Quad11%number_nonzeros + &
        qA%Quad12%number_nonzeros + &
        qA%Quad21%number_nonzeros + &
        qA%Quad22%number_nonzeros
    ENDIF

  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree

  !> Convert a dense 2nd order matrix to SpAMM matrix.
  !!
  !! @param A_dense The dense matrix.
  !!
  !! @result A The SpAMM matrix.
  function spamm_convert_dense_to_matrix_2nd_order (A_dense) result (A)

    type(spamm_matrix_2nd_order), pointer :: A
    real(spamm_kind), dimension(:, :), intent(in) :: A_dense

    A => spamm_allocate_matrix_2nd_order(size(A_dense, 1), size(A_dense, 2))
    call spamm_convert_dense_2_qutree(A_dense, A%root, 1, A%N_padded, 1, A%N_padded)

  end function spamm_convert_dense_to_matrix_2nd_order

END MODULE SpAMM_CONVERT
