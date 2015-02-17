module spamm_elementals

   use spamm_structures
   use spamm_xstructors
   use spamm_decoration

  implicit none

CONTAINS



 RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur (qA, qC)

    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC

    IF(.NOT.ASSOCIATED(qA)) RETURN

    IF(.NOT.ASSOCIATED(qC)) THEN
       CALL NewQuNode(qC, qA%i_lower, qA%j_lower, qA%i_upper, qA%j_upper)
    ENDIF

    LOG_DEBUG("q: "//to_string(qA))

    qC%Norm = qA%Norm

    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
       IF(.NOT. allocated(qC%Blok)) THEN
          ALLOCATE(qC%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
       ENDIF
       qC%Blok = qA%Blok
    ELSE
       IF(ASSOCIATED(qA%Quad11)) THEN
          !$OMP TASK UNTIED SHARED(qA,qC)
          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad11, qC%Quad11)
          !$OMP END TASK
       ENDIF
       IF(ASSOCIATED(qA%Quad12)) THEN
          !$OMP TASK UNTIED SHARED(qA,qC)
          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad12, qC%Quad12)
          !$OMP END TASK
       ENDIF
       IF(ASSOCIATED(qA%Quad21)) THEN
          !$OMP TASK UNTIED SHARED(qA,qC)
          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad21, qC%Quad21)
          !$OMP END TASK
       ENDIF
       IF(ASSOCIATED(qA%Quad22)) THEN
          !$OMP TASK UNTIED SHARED(qA,qC)
          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad22, qC%Quad22)
          !$OMP END TASK
       ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur


  SUBROUTINE SpAMM_scalar_times_tree_2d_symm_recur(beta, a)

    type(SpAMM_tree_2d_symm), pointer, intent(in)   :: A
    real(SpAMM_KIND)                                :: beta
    integer, dimension(:,:),  pointer               :: bb

    if(.not.associated(a))then
       return
    else

       bb => a%frill%bndbx 
       if(bb(1,1)-bb(0,1)==SBS)then

          a%chunk(1:SBS,1:SBS) = beta*a%chunk(1:SBS,1:SBS)
          a%frill%flops=a%frill%flops+SBS2

       else

          call SpAMM_scalar_times_tree_2d_symm_recur(beta, a%child_00)
          call SpAMM_scalar_times_tree_2d_symm_recur(beta, a%child_01)
          call SpAMM_scalar_times_tree_2d_symm_recur(beta, a%child_11)

       end if

       call SpAMM_redecorate_tree_2d_symm(a)

    end if

  END SUBROUTINE SpAMM_scalar_times_tree_2d_symm_recur


  recursive function SpAMM_absmax_tree_2d_symm (A) result(absmax)

    type(SpAMM_tree_2d_symm), pointer, intent(in)   :: A
    real(SpAMM_KIND)                                :: absmax
    integer, dimension(:,:),  pointer               :: bb

    absmax = 0
    if(.not. associated(A))then
       return
    else

       bb => A%frill%bndbx 
       if(bb(1,1)-bb(0,1)==SBS)then

          absmax = maxval(abs(a%chunk))

       else

          absmax = max(absmax, absmax_tree_2d_symm(a%child_00))
          absmax = max(absmax, absmax_tree_2d_symm(a%child_01))
          absmax = max(absmax, absmax_tree_2d_symm(a%child_11))

       endif

    endif

  end function SpAMM_absmax_tree_2d_symm




END module spamm_elemental
