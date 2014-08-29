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

  USE SpAMM_ALGEBRA
  use spamm_globals
  use spamm_management
  use spamm_types
  use spamm_utilities

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
  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree (A, qA, i_lower, i_upper, j_lower, j_upper)

    REAL(SpAMM_KIND), DIMENSION(:,:), INTENT(IN) :: A
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    INTEGER, INTENT(IN) :: i_lower, i_upper
    INTEGER, INTENT(IN) :: j_lower, j_upper

    INTEGER :: i, j
    INTEGER :: A_rows, A_cols
    integer :: convert_rows, convert_columns

    if(associated(qA)) then
      call write_log(FATAL, [ "qA should not already be associated" ])
    endif

    call write_log(2, [ "converting: "//to_string(i_lower)//" " &
      //to_string(i_upper)//" "//to_string(j_lower)//" "//to_string(j_upper) ])

    if(i_lower > size(A, 1) .or. j_lower > size(A, 2)) then
      call write_log(2, [ "outside dense matrix" ])
      return
    endif

    A_rows = i_upper-i_lower+1
    A_cols = j_upper-j_lower+1

    if(.not. associated(qA)) then
      call write_log(2, [ "allocating new node" ])
      allocate(qA)
      qA%i_lower = i_lower
      qA%i_upper = i_upper
      qA%j_lower = j_lower
      qA%j_upper = j_upper
    endif

    call write_log(2, [ "q: "//to_string(i_lower)//" "//to_string(i_upper)//" " &
      //to_string(j_lower)//" "//to_string(j_upper) ])

    IF(A_rows <= SPAMM_BLOCK_SIZE .AND. A_cols <= SPAMM_BLOCK_SIZE)THEN
      IF(A_rows < SPAMM_BLOCK_SIZE .OR. A_cols < SPAMM_BLOCK_SIZE) THEN
        call write_log(FATAL, [ "[XgpSLv6M8u5ASgg3] LOGIC ERROR IN SpAMM: padding error", &
          "A_rows = "//to_string(A_rows), &
          "A_cols = "//to_string(A_cols), &
          "SPAMM_BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE) ])
      ELSE
        call write_log(2, [ "allocating new blok" ])

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

        qA%norm = sum(matmul(qA%blok, transpose(qA%blok)), reshape( &
          (/ ((i == j, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) /), &
          (/ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE /)))

        !write(*, *) "q%blok"
        !do i = 1, SPAMM_BLOCK_SIZE
        !  write(*, "(4E10.3)") (qA%blok(i, j), j = 1, SPAMM_BLOCK_SIZE)
        !enddo

        ! Count number non-zeros.
        qA%number_nonzeros = sum(reshape( &
          (/ ((1, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) /), &
          (/ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE /)), &
          reshape( &
          (/ ((qA%blok(i, j) /= 0.0, i = 1, SPAMM_BLOCK_SIZE), j = 1, SPAMM_BLOCK_SIZE) /), &
          (/ SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE /)))

        call write_log(2, [ "non-zeros: "//to_string(qA%number_nonzeros) ])
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

      qA%number_nonzeros = 0
      qA%norm = 0

      if(associated(qA%quad11)) then
        qA%number_nonzeros = qA%number_nonzeros+qA%quad11%number_nonzeros
        qA%norm = qA%norm+qA%quad11%norm**2
      endif

      if(associated(qA%quad12)) then
        qA%number_nonzeros = qA%number_nonzeros+qA%quad12%number_nonzeros
        qA%norm = qA%norm+qA%quad12%norm**2
      endif

      if(associated(qA%quad21)) then
        qA%number_nonzeros = qA%number_nonzeros+qA%quad21%number_nonzeros
        qA%norm = qA%norm+qA%quad21%norm**2
      endif

      if(associated(qA%quad22)) then
        qA%number_nonzeros = qA%number_nonzeros+qA%quad22%number_nonzeros
        qA%norm = qA%norm+qA%quad22%norm**2
      endif

      qA%norm = sqrt(qA%norm)

    ENDIF

    call write_log(2, [ "done, going back up" ])

  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree

  !> Convert a dense 2nd order matrix to SpAMM matrix.
  !!
  !! @param A_dense The dense matrix.
  !!
  !! @result A The SpAMM matrix.
  function spamm_convert_dense_to_matrix_2nd_order (A_dense) result (A)

    type(spamm_matrix_2nd_order), pointer :: A
    real(spamm_kind), dimension(:, :), intent(in) :: A_dense

    call write_log(1, [ "converting dense matrix" ])
    A => spamm_allocate_matrix_2nd_order(size(A_dense, 1), size(A_dense, 2))
    call spamm_convert_dense_2_qutree(A_dense, A%root, 1, A%N_padded, 1, A%N_padded)

    if(associated(A%root)) then
      A%norm = A%root%norm
      A%number_nonzeros = A%root%number_nonzeros
    endif

    call write_log(1, [ "norm = "//to_string(A%norm), &
      "nnonzeros = "//to_string(A%number_nonzeros), &
      "done converting" ])

  end function spamm_convert_dense_to_matrix_2nd_order

END MODULE SpAMM_CONVERT
