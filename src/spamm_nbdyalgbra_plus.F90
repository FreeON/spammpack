!> Functions for additions.
!!
!! SpAMM generalized n-body algebras for the add
!! generalized addition (+)
!! ... [TREE-ONE-D + TREE-ONE-D] ... [TREE-ONE-D + TREE-ONE-D] ...
!!   SpAMM_tree_1d_plus_tree_1d
module spamm_nbdyalgbra_plus

  use spamm_decoration
  use spamm_elementals
  use spamm_nbdyalgbra_times
  use spamm_structures
  use spamm_xstructors

  implicit none

contains

  !++NBODYADDTN:     c => inplace*c + beta*b
  function SpAMM_tree_1d_plus_tree_1d (C, beta, B, inplace) result(D)

    type(SpAMM_tree_1d), pointer, intent(IN)    :: B
    type(SpAMM_tree_1d), pointer, intent(INOUT) :: C
    real(SPAMM_KIND),             intent(IN)    :: Beta
    real(SPAMM_KIND),    optional,intent(in)    :: inplace
    type(SpAMM_tree_1d), pointer                :: d
    real(SPAMM_KIND)                            :: Alpha
    !
    D=>C
    !
    if(.not. associated(B))return
    if(beta==0)   return
    !
    if(present(inplace))then
       alpha=inplace
    else
       alpha=SpAMM_one
    end if

    call SpAMM_tree_1d_plus_tree_1d_inplace_recur(d, B, alpha, beta )

  end function SpAMM_tree_1d_plus_tree_1d

  ! for tree_1d, A = alpha*A + beta*B
  recursive subroutine SpAMM_tree_1d_plus_tree_1d_inplace_recur(a, b, alpha, beta)

    type(SpAMM_tree_1d), pointer                :: A
    type(SpAMM_tree_1d), pointer, intent(IN)    :: B
    real(SPAMM_KIND),             intent(IN)    :: alpha, beta
    logical                                     :: TA, TB

    TA=associated(A)
    TB=associated(B)

    if(.not.TA.and..not.TB)then

       return

    elseif(TA .and. .not.TB) then

       ! a=alpha*a
       call SpAMM_scalar_times_tree_1d_recur(alpha, a)

    elseif(TA.and.TB)then

       if(b%frill%leaf)then
          ! A = alpha*A + beta*B

          a%frill%needs_initialization=.false.
          a%chunk(1:SBS) = alpha*a%chunk(1:SBS) + beta*b%chunk(1:SBS)
          a%frill%flops = a%frill%flops + 3*SBS

       else ! recursively descend ...

          call SpAMM_tree_1d_plus_tree_1d_inplace_recur(spamm_construct_child_1d(a, 0), & !0>
               b%child_0, alpha, beta)

          call SpAMM_tree_1d_plus_tree_1d_inplace_recur(spamm_construct_child_1d(a, 0), & !1>
               b%child_1, alpha, beta)
       end if

    end if

    call SpAMM_redecorate_tree_1d(a)

  end subroutine SpAMM_tree_1d_plus_tree_1d_inplace_recur


  !!
  !! ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ...
  !!


  !> Add a scaled identity matrix (in place).
  !!
  !! @f$ d \leftarrow a + \alpha \,\, \mathrm{Id} @f$.
  !!
  !! @param alpha The scaling factor.
  !! @param a The matrix.
  !! @return The matrix sum.
  function spamm_scalar_plus_tree_2d_symm(alpha, a) result(d)

    real(SPAMM_KIND), intent(in) :: alpha
    type(spamm_tree_2d_symm), pointer :: a
    type(spamm_tree_2d_symm), pointer :: d

    if(.not. associated(a)) then
       d => null()
       return
    else
       d => a
       call spamm_scalar_plus_tree_2d_symm_recur(alpha, d)
    end if

  end function spamm_scalar_plus_tree_2d_symm

  recursive subroutine spamm_scalar_plus_tree_2d_symm_recur(alpha, a)

    real(SPAMM_KIND), intent(in) :: alpha
    type(spamm_tree_2d_symm), pointer :: a

    integer, dimension(2) :: lo,hi
    integer :: i

    if(.not. associated(a)) return
    if(a%frill%leaf) then
       a%frill%needs_initialization = .false.
       lo = a%frill%bndbx(0, :)
       hi = a%frill%bndbx(1, :)
       do i = 1, hi(1)-lo(1)+1
          a%chunk(i, i) = a%chunk(i, i)+alpha
       end do
    else
       ! child along [00]:
       call spamm_scalar_plus_tree_2d_symm_recur(alpha, spamm_construct_tree_2d_symm_00(a))
       ! child along [11]:
       call spamm_scalar_plus_tree_2d_symm_recur(alpha, spamm_construct_tree_2d_symm_11(a))
    end if
    call spamm_redecorate_tree_2d_symm(a)

  end subroutine spamm_scalar_plus_tree_2d_symm_recur

  function spamm_tree_2d_symm_plus_tree_2d_symm(a, b, alpha, beta, c) result(d)

    type(spamm_tree_2d_symm), pointer, intent(inout)           :: a, b
    type(spamm_tree_2d_symm), pointer, intent(inout), optional :: c
    real(spamm_kind), intent(in), optional                     :: alpha, beta

    type(spamm_tree_2d_symm), pointer                          :: d

    real(spamm_kind)                                           :: local_alpha,local_beta

    D=>NULL()

    if(.not. associated(A))return
    if(.not. associated(B))return

    if(present(Alpha))then; Local_Alpha=Alpha; else; Local_Alpha=SpAMM_One;end if
       if(present(Beta ))then; Local_Beta =Beta;  else; Local_Beta=SpAMM_One; end if

          if(present(C))then ! we are going for an in place plus with an existing C:

             if(associated(B,C))then  ! if passed in C is B, then in place accumulation on B ...

                ! B => alpha*A + beta*B
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(B, A, Local_beta, Local_alpha)
                D=>B

             elseif(associated(A,C))then ! if passed in C is A, then in place accumulation on A ...

                ! A => alpha*A + beta*B ...
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(A, B, Local_alpha, Local_beta)
                D=>A

             else  ! C is passed in as a seperate channel for accumulation ...

                ! C => C + alpha*A + beta*B
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(C, A, B, Local_alpha, Local_beta)
                D=>C

             end if

          else  ! We need a clean tree_2d_symm at this point ...

             D => SpAMM_new_top_tree_2d_symm(A%frill%NDimn)
             ! D => D + alpha*A + beta*B
             call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(C, A, B, Local_alpha, Local_beta)
             D=>C

          end if

          call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(D, A, B, alpha, beta)

        end function spamm_tree_2d_symm_plus_tree_2d_symm

        ! for tree_2d_symm, A = A + alpha*A + beta*B
        recursive subroutine SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(a, b, alpha, beta)

          type(SpAMM_tree_2d_symm), pointer                :: A
          type(SpAMM_tree_2d_symm), pointer, intent(IN)    :: B
          real(SPAMM_KIND),                  intent(IN)    :: alpha, beta
          logical                                          :: TA, TB

          TA=associated(A)
          TB=associated(B)

          if(.not.TA.and.TB)then

             ! a = b (zero thresholding)
             call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur(a, b, 0.0_SPAMM_KIND)

             ! a = beta*b = beta*a
             !       CALL SpAMM_scalar_times_tree_2d_symm_recur(beta, a)

          elseif(TA .and. .not.TB) then

             ! a=alpha*a
             !       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a)

          elseif(TA.and.TB)then

             if(b%frill%leaf)then

                ! A = alpha*A + beta*B
                a%frill%needs_initialization=.false.
                a%chunk(1:SBS,1:SBS) = alpha*a%chunk(1:SBS,1:SBS) + beta*b%chunk(1:SBS,1:SBS)
                a%frill%flops = a%frill%flops + 3*SBS2

             else

                ! recursively decend, possibly building out A if nessesary ...
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_00(a), & !00>
                     b%child_00, alpha, beta)
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_01(a), & !01>
                     b%child_01, alpha, beta)
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_10(a), & !10>
                     b%child_10, alpha, beta)
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(SpAMM_construct_tree_2d_symm_11(a), & !11>
                     b%child_11, alpha, beta)
             end if

             ! enrich the resultnt
             call SpAMM_redecorate_tree_2d_symm(a)

          end if

        end subroutine SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur

        ! for tree_2d_symm: C = C + alpha*A+beta*B
        recursive subroutine SpAMM_tree_2d_symm_plus_tree_2d_symm_recur(C, A, B, alpha, beta)

          type(SpAMM_tree_2d_symm), pointer, intent(IN)    :: A,B
          type(SpAMM_tree_2d_symm), pointer                :: C
          real(SPAMM_KIND)                                 :: alpha, beta
          logical                                          :: TA, TB

          TA=associated(A)
          TB=associated(B)

          if(TA .and. .not.TB) then

             ! C = C + alpha*A
             call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(c, a, alpha, SpAMM_One)

          elseif(.not.TA.and.TB)then

             ! C = C + beta*B
             call SpAMM_tree_2d_symm_plus_tree_2d_symm_inplace_recur(c, b, SpAMM_One, beta)

          elseif(TA.and.TB)then

             ! leaf situation ...
             if(c%frill%leaf)then

                ! c = c + alpha*a + beta*b
                c%frill%needs_initialization=.false.
                c%chunk(1:sbs,1:sbs)=c%chunk(1:sbs,1:sbs)+alpha*a%chunk(1:sbs,1:sbs)+beta*b%chunk(1:sbs,1:sbs)
                c%frill%flops=c%frill%flops+3*sbs2

             else
                ! recursively decend, popping new children as needed ...
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_00(c), & !00>
                     a%child_00, b%child_00, alpha, beta)
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_01(c), & !01>
                     a%child_01, b%child_01, alpha, beta)
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_10(c), & !10>
                     a%child_10, b%child_10, alpha, beta)
                call SpAMM_tree_2d_symm_plus_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_11(c), & !11>
                     a%child_11, b%child_11, alpha, beta)
             end if

             call SpAMM_redecorate_tree_2d_symm(c)

          end if

        end subroutine SpAMM_tree_2d_symm_plus_tree_2d_symm_recur

      end module spamm_nbdyalgbra_plus
