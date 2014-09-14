!> Copyright (c) 2014, Los Alamos National Laboratory
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

!> Multiply a vector with a scalar: @f$ V \leftarrow \alpha V @f$.
!!
!! @param V The vector
!! @param alpha The scalar @f$ \alpha @f$
subroutine spamm_multiply_order_1_x_scalar (V, alpha)

  type(spamm_matrix_order_1), pointer, intent(inout) :: V
  real(SPAMM_KIND), intent(in) :: alpha

  call spamm_multiply_bitree_x_scalar(V%root, alpha)

end subroutine spamm_multiply_order_1_x_scalar

!> Scalar multiply: \f$ A \leftarrow \alpha A \f$.
!!
!! @param bA Pointer to vector A.
!! @param alpha Scalar @f$ \alpha @f$.
RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar(bA, alpha)

  TYPE(BiTree), POINTER :: bA
  REAL(SpAMM_KIND) :: alpha

  IF(.NOT.ASSOCIATED(bA))RETURN

  !$OMP TASK SHARED(bA)
  CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA, alpha)
  !$OMP END TASK
  !$OMP TASKWAIT

END SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar

!> Allocate a 1st order SpAMM matrix.
!!
!! @param N The number of columns.
!! @param V The new vector.
subroutine spamm_allocate_matrix_order_1 (N, V)

  use spamm_globals

  integer, intent(in) :: N
  type(spamm_matrix_order_1), pointer, intent(out) :: V

  V => null()
  allocate(V)

  V%N = N
  V%depth = 0
  V%N_padded = SPAMM_BLOCK_SIZE

  ! This should be pretty efficient for reasonable matrix sizes and is presumably faster than some logarithm calculation since it
  ! only involves an increment and a bit shift.
  do while(V%N_padded < N)
    !LOG_DEBUG("depth = "//to_string(V%depth)//", N_padded = "//to_string(V%N_padded))
    V%depth = V%depth+1
    V%N_padded = 2*V%N_padded
  enddo

  !LOG_INFO("allocated "//to_string(N)//" vector")
  !LOG_INFO("  BLOCK_SIZE = "//to_string(SPAMM_BLOCK_SIZE))
  !LOG_INFO("  N_padded   = "//to_string(V%N_padded))
  !LOG_INFO("  depth      = "//to_string(V%depth))

end subroutine spamm_allocate_matrix_order_1
