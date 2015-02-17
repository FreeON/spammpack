module spamm_xstructors

  use spamm_structures
  use spamm_decoration

  implicit none

contains

  ! SpAMM xstructors ________________ SXTRS _________________
  ! constructors and destructors 
  !
  ! instantiates the tree_2d_symm structure 

  function SpAMM_new_top_tree_2d_symm ( M, N ) result (tree)
    !
    integer,                intent(in) :: M, N
    integer                            :: M_pad, N_pad, depth
    type(SpAMM_tree_2d_symm), pointer  :: tree

    ! instantiate the root node.  this is the tree top ... 
    allocate(tree)

    ! here are padded dimensions ...
    do depth=0,64
       M_pad=SPAMM_BLOCK_SIZE*2**depth
       if(M_pad>M)exit
    enddo
    !
    do depth=0,64
       N_pad=SPAMM_BLOCK_SIZE*2**depth
       if(N_pad>N)exit
    enddo

    ! the [i-j] native dimensions ...
    tree%frill%ndimn        = (/ M, N /)    

    ! and the bounding box
    allocate( tree%frill%bndbx(0:1,1:2) )

    ! the padded tiles
    tree%frill%bndbx(0:1,1) = (/ 0, M_pad /)  ! [i-lo,i-hi]
    tree%frill%bndbx(0:1,2) = (/ 0, N_pad /)  ! [j-lo,j-hi]

  end function SpAMM_new_top_tree_2d_symm
  
  ! structor for the lo-lo [00] channel ...
  function SpAMM_construct_tree_2d_00(tree) result(ch00)  

    type(SpAMM_tree_2d_symm), intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer       :: ch00
    integer, dimension(:,:),  pointer       :: bb, bb00
    integer, dimension(1:2)                 :: mid,wid

    ch00=>tree%child_00
    ! pre-existing?  ok, later
    if(associated(ch00))return

    ! otherwise, instantiate ...
    allocate(ch00)

    ! global... do we really need this? 
    ch00%frill%ndimn = tree%frill%ndimn

    ! local boxes ... 
    bb=>tree%frill%bndbx
    bb00=>ch00%frill%bndbx
    if(.not.associated(bb00))allocate(bb00(0:1,1:2))

    ! the 00 split ...
    wid(:)=bb(1,:)-bb(0,:) 
    mid(:)=bb(0,:)+wid(:)/2 
    bb00(:,1)=(/ bb(0,1) , mid(1) /)
    bb00(:,2)=(/ bb(0,2) , mid(2) /)

  end function SpAMM_construct_tree_2d_00

  ! structor for the lo-hi [01] channel ...
  function SpAMM_construct_tree_2d_01(tree) result(ch01)  

    type(SpAMM_tree_2d_symm), intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer       :: ch01
    integer, dimension(:,:) , pointer       :: bb, bb01
    integer, dimension(1:2)                 :: mid,wid

    ch01=>tree%child_01
    ! pre-existing?  ok, later
    if(associated(ch01))return

    ! otherwise, instantiate ...
    allocate(ch01)

    ! global... do we really need this? 
    ch01%frill%ndimn = tree%frill%ndimn

    ! local boxes ... 
    bb=>tree%frill%bndbx
    bb01=>ch01%frill%bndbx
    if(.not.associated(bb01))allocate(bb01(0:1,1:2))

    ! the 01 split ...
    wid(:)=bb(1,:)-bb(0,:) 
    mid(:)=bb(0,:)+wid(:)/2 
    bb01(:,1)=(/bb(0,1) , mid(1) /)
    bb01(:,2)=(/mid(2)+1, bb(1,2)/)

  end function SpAMM_construct_tree_2d_01

  ! structor for the hi-hi [11] channel ...
  function SpAMM_construct_tree_2d_11(tree) result(ch11)  

    type(SpAMM_tree_2d_symm), intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer       :: ch11
    integer, dimension(:,:),  pointer       :: bb, bb11
    integer, dimension(1:2)                 :: mid,wid

    ch11=>tree%child_11
    ! pre-existing?  ok, later
    if(associated(ch11))return

    ! otherwise, instantiate ...
    allocate(ch11)

    ! global... do we really need this? 
    ch11%frill%ndimn = tree%frill%ndimn

    ! local boxes ... 
    bb=>tree%frill%bndbx
    bb11=>ch11%frill%bndbx
    if(.not.associated(bb11))allocate(bb11(0:1,1:2))

    ! the 11 split ...
    wid(:)=bb(1,:)-bb(0,:) 
    mid(:)=bb(0,:)+wid(:)/2 
    bb11(0,:)=mid(:)+1     ! hi,hi [mid+1, hi]
    bb11(1,:)=bb(1,:)      ! hi,hi [mid+1, hi]

  end function SpAMM_construct_tree_2d_11

  ! Structor to recursively destroy a tree_2d_symm ... 
  recursive subroutine  SpAMM_destruct_tree_2d_symm_recur (self)
    !
    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: self
    type(SpAMM_tree_2d_symm), pointer                :: ch ! sub-tree pointer 
 
    ! check for self-non-association (eg. at leaf pntr) ...
    if(.not.associated(self))return
    
    ch=>self%child_00                            ! take the [00] channel
    call SpAMM_destruct_tree_2d_symm_recur (ch)  ! recur to the botom:    
    call SpAMM_destruct_tree_2d_symm_node  (ch)  ! kill backwards up the tree
    ch=>self%child_01                            ! take the [01] channel
    call SpAMM_destruct_tree_2d_symm_recur (ch)  ! recur to the botom:    
    call SpAMM_destruct_tree_2d_symm_node  (ch)  ! kill backwards up the tree 
    ch=>self%child_11                            ! take the [11] channel
    call SpAMM_destruct_tree_2d_symm_recur (ch)  ! recur to the botom:    
    call SpAMM_destruct_tree_2d_symm_node  (ch)  ! kill backwards up the tree 

  end subroutine SpAMM_destruct_tree_2d_symm_recur ! ... and we're out ... 

  subroutine  SpAMM_destruct_tree_2d_symm_node (self)

    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: self

    ! kill the adornment, mort le accoutrement ...    
    if(allocated(self%chunk))deallocate(self%chunk) 
    deallocate(self%frill%bndbx)   ! fru-fru 
    deallocate(self)               ! done 
    nullify(self)                  ! bye-bye

  end subroutine SpAMM_destruct_tree_2d_symm_node

  function SpAMM_new_identity_tree_2d_symm ( M, N, Alpha_O ) result (tree)
    !
    integer,        intent(in)  :: M, N
    REAL(SpAMM_KIND), OPTIONAL  :: Alpha_O
    REAL(SpAMM_KIND)            :: Alpha
    integer                     :: depth
    type(SpAMM_tree_2d_symm), pointer :: tree

    tree => SpAMM_new_top_tree_2d_symm ( M, N )

    IF(PRESENT(alpha_O))THEN
       alpha=Alpha_O
    ELSE
       alpha=SpAMM_One
    ENDIF    

    ! push alpha onto the trace ...
    depth=0
    CALL SpAMM_new_identity_tree_2d_symm_recur (tree, alpha, depth)

  end function SpAMM_new_identity_tree_2d_symm

  ! putting alpha down, onto the trace of this tree_2d_symm ...
  recursive subroutine SpAMM_new_identity_tree_2d_symm_recur (tree, alpha, depth)

    type(SpAMM_tree_2d_symm)                :: tree
    real(SpAMM_KIND),        intent(in)     :: alpha
    integer,                 intent(in)     :: depth 
    INTEGER, DIMENSION(:,:), pointer        :: bb
    integer                                 :: i 

    bb => tree%frill%bndbx 

    IF( bb(1,1)-bb(0,1) == SBS )THEN  ! Leaf condition ? 

       ! here is a 2d chunk (ch) ... 
       allocate( tree%chunk( SBS, SBS ) ) 

       ! set its trace with scalar alpha ...
       tree%chunk=SpAMM_Zero
       do i=1, SBS
          tree%chunk(i,i)=alpha 
       enddo

    ELSEIF(depth>16)THEN ! build the identity tree down

       ! child along [00]: [lo,mid]x[lo,mid] ... 
       CALL SpAMM_new_identity_tree_2d_symm_recur( SpAMM_construct_tree_2d_00(tree) , alpha, depth+1 )
       ! nothing off diagonal 
       tree%child_01=>Null() 
       ! child along [11]: [mid+1,hi]x[mid+1,hi] ... 
       CALL SpAMM_new_identity_tree_2d_symm_recur( SpAMM_construct_tree_2d_11(tree) , alpha, depth+1 )

    ELSE; STOP ' depth 16 exceeded in SpAMM_init_ident_tree_2d_symm_recur';  
    ENDIF
    
    ! merge & regarnish back up the tree ...
    CALL SpAMM_redecorate_tree_2d_symm(tree)
    !
  end subroutine SpAMM_new_identity_tree_2d_symm_recur

end module spamm_xstructors

!!$ RECURSIVE SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur (qA, qC)
!!$
!!$    TYPE(QuTree), POINTER, INTENT(IN)    :: qA
!!$    TYPE(QuTree), POINTER, INTENT(INOUT) :: qC
!!$
!!$    IF(.NOT.ASSOCIATED(qA)) RETURN
!!$
!!$    IF(.NOT.ASSOCIATED(qC)) THEN
!!$       CALL NewQuNode(qC, qA%i_lower, qA%j_lower, qA%i_upper, qA%j_upper)
!!$    ENDIF
!!$
!!$    LOG_DEBUG("q: "//to_string(qA))
!!$
!!$    qC%Norm = qA%Norm
!!$
!!$    IF(qA%i_upper-qA%i_lower+1 == SPAMM_BLOCK_SIZE) then
!!$       IF(.NOT. allocated(qC%Blok)) THEN
!!$          ALLOCATE(qC%Blok(SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE))
!!$       ENDIF
!!$       qC%Blok = qA%Blok
!!$    ELSE
!!$       IF(ASSOCIATED(qA%Quad11)) THEN
!!$          !$OMP TASK UNTIED SHARED(qA,qC)
!!$          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad11, qC%Quad11)
!!$          !$OMP END TASK
!!$       ENDIF
!!$       IF(ASSOCIATED(qA%Quad12)) THEN
!!$          !$OMP TASK UNTIED SHARED(qA,qC)
!!$          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad12, qC%Quad12)
!!$          !$OMP END TASK
!!$       ENDIF
!!$       IF(ASSOCIATED(qA%Quad21)) THEN
!!$          !$OMP TASK UNTIED SHARED(qA,qC)
!!$          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad21, qC%Quad21)
!!$          !$OMP END TASK
!!$       ENDIF
!!$       IF(ASSOCIATED(qA%Quad22)) THEN
!!$          !$OMP TASK UNTIED SHARED(qA,qC)
!!$          CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qA%Quad22, qC%Quad22)
!!$          !$OMP END TASK
!!$       ENDIF
!!$    ENDIF
!!$
!!$  END SUBROUTINE SpAMM_Copy_QuTree_2_QuTree_Recur
