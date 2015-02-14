
  ! SpAMM STuCtuReS ________________ SSTCRS _________________
  ! SpAMM constructers and destructures -- structures


  function SpAMM_new_tree_top_2d_symm ( M, N, Alpha ) result (tree)
    !
    integer,        intent(in)  :: M, N
    REAL(SpAMM_KIND), OPTIONAL  :: Alpha_O
    REAL(SpAMM_KIND)            :: Alpha
    integer                     :: M_pad, N_pad
    type(tree_2d_symm), pointer :: Tree

    ! the root node (top) of the tree ... 
    allocate(tree)

    ! here are the padded dimensions ...
    do depth=0,64
       M_pad=SPAMM_BLOCK_SIZE*2**depth
       if(M_pad>M)exit
    enddo
    !
    do depth=0,64
       N_pad=SPAMM_BLOCK_SIZE*2**depth
       if(N_pad>N)exit
    enddo

    ! the tree dimensions in frills ...
    tree%frill%ndimn = (/ M, N /)
    tree%frill%bndbx = ( (/     0,     0 /) , &
                         (/ M_pad, N_pad /) )

    IF(PRESENT(alpha_O))THEN
       alpha=Alpha_O
    ELSE
       alpha=SpAMM_Zero
    ENDIF

    tree%child_01=>Null()

    depth=0
    CALL SpAMM_init_ident_2d_symmetric (tree%child_00, alpha, depth)

    depth=0
    CALL SpAMM_init_ident_2d_symmetric (tree%child_11, alpha, depth)

  end function SpAMM_init_ident_tree_2d_symm

  recursive subroutine SpAMM_init_ident_tree_2d_symm_recur (tree, alpha, depth)

    use spamm_bisect
    use spamm_globals

    type(tree_2d_symmetric), intent(inout) :: tree
    real(SpAMM_KIND),        intent(in)    :: alpha
    integer,                 intent(in)    :: depth 
    ! 
    INTEGER, DIMENSION(2,2), pointer       :: bb, bb00, bb11

    !
    bb => tree%frill%bndbx

    ! w = [1]-[0]
    mwid=bb(1,1)-bb(0,1)
    
    IF( mwid == SPAMM_BLOCK_SIZE )THEN  ! Leaf condition ? 

       ! here is a 2d chunk (ch) ... 
       allocate( tree%chunk( SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE ) );  ch => tree%chunk

       ! set its trace with scalar alpha ...
       ch=SpAMM_Zero
       do i=1, SPAMM_BLOCK_SIZE; ch(i:i)=alpha; enddo

    ELSEIF(depth>16)THEN ! build the identity tree down

       ! w = [1]-[0]
       nwid=bb(1,2)-bb(0,2) ! =mwid 

       ! mid=lo+wid/2 in [ixj] ...
       mmid=bb(0,1)+mwid/2 ! i
       nmid=bb(0,2)+nwid/2 ! j

       ! childrens (building down right hyar) ...
       allocate(tree%child_00); allocate(tree%child_00%bndbx); bb00=>tree%child_00%frill%bndbx
       allocate(tree%child_11); allocate(tree%child_11%bndbx); bb11=>tree%child_11%frill%bndbx

       ! [00]: [lo,mid]x[lo,mid] ... 
       bb00(0,1)=bb(0,1) ! [lo,i]
       bb00(0,2)=bb(0,2) ! [lo,j]
       bb00(1,1)=mmid    ! [hi,i]
       bb00(1,2)=nmid    ! [hi,j]
       CALL SpAMM_init_ident_tree_2d_symm_recur(tree%child_00, alpha, depth+1)

       ! [11]: [mid+1,hi]x[mid+1,hi] ... 
       bb11(0,1)=mmid+1  ! [lo,i]
       bb11(0,2)=nmid+1  ! [lo,j]
       bb11(1,1)=bb(1,1) ! [hi,i]
       bb11(1,2)=bb(1,2) ! [hi,j]
       CALL SpAMM_init_ident_tree_2d_symm_recur(tree%child_11, alpha, depth+1)
       
       ! traversal enrichment ...
       !       CALL decorate(tree)
       !       
    ELSE; STOP ' depth exceeded in SpAMM_init_ident_tree_2d_symm_recur' 
    ENDIF
    !
  end subroutine SpAMM_init_ident_tree_2d_symm_recur

  !> The destructor.
  recursive subroutine  SpAMM_delete_tree_2d_symm_recur (self)
    !
    type(tree_2d_symm), pointer, intent(inout) :: self
    type(tree_2d_symm), pointer                :: ch ! sub-tree pointer 
 
    ! check for self-non-association (at leaf) ...
    if(.not.associated(self))return

    ! otherwise recur to the botom:    
    ch=>self%child_00
    call SpAMM_delete_tree_2d_symm_recur (ch)
    ! backwards kill of 00 and accoutrement ...
    if(allocated(ch%chunk))deallocate(ch%chunk)
    deallocate(ch%frill%bndbx)   ! fru-fru
    deallocate(ch)               ! done 
    nullify(ch)                  ! bye-bye

    ch=>self%child_01
    call SpAMM_delete_tree_2d_symm_recur (ch)
    ! backwards kill of 01 and accoutrement ...    
    if(allocated(ch%chunk))deallocate(ch%chunk)
    deallocate(ch%frill%bndbx)   ! fru-fru
    deallocate(ch)               ! done 
    nullify(ch)                  ! bye-bye

    ch=>self%child_11
    call SpAMM_delete_tree_2d_symm_recur (ch)
    ! backwards kill of 11 and accoutrement ...    
    if(allocated(ch%chunk))deallocate(ch%chunk)
    deallocate(ch%frill%bndbx)   ! fru-fru
    deallocate(ch)               ! done 
    nullify(ch)                  ! bye-bye

    ! step back up the tree 
  end subroutine SpAMM_delete_tree_2d_symm_recur ! ... and we're out ... 
  !  

