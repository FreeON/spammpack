!> Management functions for quadtree objects.
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
!! @author Nicolas Bock nicolasbock@freeon.org
module spamm_management_1d

  use spamm_types
  use spamm_real_precision

  implicit none

  private

contains

  !> Copy a vector to another vector:
  !!
  !! \f$ B \leftarrow A \f$.
  !!
  !! @param A The vector A.
  !! @param B The vector B.
  !> recursive subroutine copy_matrix_1d_to_matrix_1d (A, B)

  !>   type(spamm_matrix_1d), pointer, intent(in) :: A
  !>   type(spamm_matrix_1d), pointer, intent(inout) :: B

  !>   if(.not. associated(A)) then
  !>      return
  !>   endif

  !>   if(associated(B)) then
  !>      call delete(B)
  !>      B => null()
  !>   endif

  !>   LOG_DEBUG("copying vector")

  !>   call new(A%N, B)

  !>   if(.not. associated(bC)) then
  !>      call NewBiNode(bC, bA%i_lower, bA%i_upper)
  !>   endif

  !>   bC%Norm=bA%Norm
  !>   if(bA%i_upper-bA%i_lower+1 == SPAMM_BLOCK_SIZE .and. allocated(bA%vect)) then
  !>      if(.not. allocated(bC%Vect)) then
  !>         allocate(bC%Vect(SPAMM_BLOCK_SIZE))
  !>      endif
  !>      LOG_DEBUG("copying vect")
  !>      bC%Vect = bA%Vect
  !>   ELSE
  !>      IF(ASSOCIATED(bA%sect1))THEN
  !>         !$OMP TASK UNTIED SHARED(bA,bC) &
  !>         !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
  !>         CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA%sect1,bC%sect1,Depth+1)
  !>         !$OMP END TASK
  !>      ELSEIF(ASSOCIATED(bC%sect1))THEN
  !>         !$OMP TASK UNTIED SHARED(bC) &
  !>         !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
  !>         CALL SpAMM_Delete_BiTree_Recur(bC%sect1)
  !>         !$OMP END TASK
  !>         !$OMP TASKWAIT
  !>         DEALLOCATE(bC%sect1)
  !>      ENDIF
  !>      IF(ASSOCIATED(bA%sect2))THEN
  !>         !$OMP TASK UNTIED SHARED(bA,bC) &
  !>         !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
  !>         CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA%sect2,bC%sect2,Depth+1)
  !>         !$OMP END TASK
  !>      ELSEIF(ASSOCIATED(bC%sect2))THEN
  !>         !$OMP TASK UNTIED SHARED(bC) &
  !>         !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
  !>         CALL SpAMM_Delete_BiTree_Recur(bC%sect2)
  !>         !$OMP END TASK
  !>         !$OMP TASKWAIT
  !>         DEALLOCATE(bC%sect2)
  !>      ENDIF
  !>   ENDIF

  !> end subroutine copy_matrix_1d_to_matrix_1d

end module spamm_management_1d
