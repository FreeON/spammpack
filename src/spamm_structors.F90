
module spamm_structers

  use spamm_globals
  use spamm_types
  implicit none

contains

  ! SpAMM structer ________________ STCTRS _________________
  ! constructers and destructers 
  !
  ! instantiates the tree_2d_symm structure 

  function SpAMM_new_top_tree_2d_symm ( M, N ) result (tree)
    !
    integer,        intent(in)  :: M, N
    integer                     :: M_pad, N_pad
    type(tree_2d_symm), pointer :: tree

    ! instantiate the root node, tree top ... 
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

    ! the [i-j] tree dimensional frill ...
    tree%frill%ndimn = (/ M, N /)    
    tree%frill%bndbx(0:1,1) = (/ 0, M_pad /)  ! [i-lo,i-hi]
    tree%frill%bndbx(0:1,2) = (/ 0, N_pad /)  ! [j-lo,j-hi]

  end function SpAMM_new_top_tree_2d_symm

  function SpAMM_set_new_identity_tree_2d_symm ( M, N, Alpha ) result (tree)
    !
    integer,        intent(in)  :: M, N
    REAL(SpAMM_KIND), OPTIONAL  :: Alpha_O
    REAL(SpAMM_KIND)            :: Alpha
    integer                     :: M_pad, N_pad
    type(tree_2d_symm), pointer :: tree

    tree = SpAMM_new_top_tree_2d_symm ( M, N )

    IF(PRESENT(alpha_O))THEN
       alpha=Alpha_O
    ELSE
       alpha=SpAMM_One
    ENDIF    

    ! push alpha on 00 down to the trace 
    depth=0
    CALL SpAMM_init_ident_tree_2d_symm_recur (tree%child_00, alpha, depth)
    ! kill the diagonal 
    !    CALL SpAMM_delete_tree_2d_symm_recur     (tree%child_01)
    ! push alpha on 11 down to the trace 
    depth=0
    CALL SpAMM_init_ident_tree_2d_symm_recur (tree%child_11, alpha, depth)

  end function SpAMM_set_new_identity_tree_2d_symm

  ! putting alpha down, onto the trace of this tree_2d_symm ...
  recursive subroutine SpAMM_init_ident_tree_2d_symm_recur (tree, alpha, depth)

    use spamm_bisect
    use spamm_globals

    type(tree_2d_symmetric), intent(inout) :: tree
    real(SpAMM_KIND),        intent(in)    :: alpha
    integer,                 intent(in)    :: depth 
    INTEGER, DIMENSION(2,2), pointer       :: bb, bb00, bb11

    bb => tree%frill%bndbx 
    wid=bb(1,1)-bb(0,1) ! w = [i]-[o]   

    IF( wid == SPAMM_BLOCK_SIZE )THEN  ! Leaf condition ? 

       ! here is a 2d chunk (ch) ... 
       allocate( tree%chunk( SPAMM_BLOCK_SIZE, SPAMM_BLOCK_SIZE ) ) 
       ! set its trace with scalar alpha ...
       tree%chunk=SpAMM_Zero
       do i=1, SPAMM_BLOCK_SIZE; tree%chunk(i:i)=alpha; enddo

    ELSEIF(depth>16)THEN ! build the identity tree down

       ! child along [00]: [lo,mid]x[lo,mid] ... 
       CALL SpAMM_init_ident_tree_2d_symm_recur( SpAMM_tree_2d_construct_00(tree) , alpha, depth+1 )
       tree%child_01=>Null() ! nothing off diagonal 
       ! child along [11]: [mid+1,hi]x[mid+1,hi] ... 
       CALL SpAMM_init_ident_tree_2d_symm_recur( SpAMM_tree_2d_construct_11(tree) , alpha, depth+1 )

    ELSE; STOP ' depth 16 exceeded in SpAMM_init_ident_tree_2d_symm_recur';  
    ENDIF
    
    ! merge & regarnish back up the tree ...
    CALL SpAMM_redecorate_2d_symm(tree)
    !
  end subroutine SpAMM_init_ident_tree_2d_symm_recur
  
  ! structer for the lo-lo [00] channel ...
  function SpAMM_tree_2d_construct_00(tree) result(ch00)  

    type(SpAMM_tree_2d_symm), intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer       :: ch00
    integer, dimension(0:1,1:2), pointer    :: bb, bb00
    integer, dimension(1:2)                 :: mid,wid

    ! instantiate ...
    ch00=>tree%child_00
    if(.not.associated(ch00))then
       allocate(ch00)
       allocate(ch00%frill%bndbx)
    endif

    ! global... do we really need this? 
    ch%frill%ndimn = tree%frill%ndimn

    ! local boxes ... 
    bb=>tree%frill%bndbx
    bb00=>ch00%frill%bndbx

    ! the split ...
    wid(:)=bb(1,:)-bb(0,:) 
    mid(:)=bb(0,:)+wid(:)/2 
    bb00(:,1)=(/ bb(0,1) , mid(1) /)
    bb00(:,2)=(/ bb(0,2) , mid(2) /)

  end function SpAMM_tree_2d_construct_00

  ! structer for the lo-hi [01] channel ...
  function SpAMM_tree_2d_construct_01(tree) result(ch01)  

    type(SpAMM_tree_2d_symm), intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer       :: ch01
    integer, dimension(0:1,1:2), pointer    :: bb, bb01
    integer, dimension(1:2)                 :: mid,wid

    ! instantiate ...
    ch01=>tree%child_01
    if(.not.associated(ch01))then
       allocate(ch01)
       allocate(ch01%frill%bndbx)
    endif

    ! global... do we really need this? 
    ch%frill%ndimn = tree%frill%ndimn

    ! local boxes ... 
    bb=>tree%frill%bndbx
    bb01=>ch01%frill%bndbx

    ! the split ...
    wid(:)=bb(1,:)-bb(0,:) 
    mid(:)=bb(0,:)+wid(:)/2 
    bb01(:,1)=(/bb(0,1) , mid(1) /)
    bb01(:,2)=(/mid(2)+1, bb(1,2)/)

  end function SpAMM_tree_2d_construct_01

  ! structer for the hi-hi [11] channel ...
  function SpAMM_tree_2d_construct_11(tree) result(ch11)  

    type(SpAMM_tree_2d_symm), intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer       :: ch11
    integer, dimension(0:1,1:2), pointer    :: bb, bb11
    integer, dimension(1:2)                 :: mid,wid

    ! instantiate ...    
    ch11=>tree%child_11
    if(.not.associated(ch11))then
       allocate(ch11)
       allocate(ch11%frill%bndbx)
    endif

    ! global... do we really need this? 
    ch%frill%ndimn = tree%frill%ndimn

    ! local boxes ... 
    bb=>tree%frill%bndbx
    bb11=>ch11%frill%bndbx

    ! the split ...
    wid(:)=bb(1,:)-bb(0,:) 
    mid(:)=bb(0,:)+wid(:)/2 
    bb11(0,:)=mid(:)+1     ! hi,hi [mid+1, hi]
    bb11(1,:)=bb(1,:)      ! hi,hi [mid+1, hi]

  end function SpAMM_tree_2d_construct_11

  ! Structer to void tree_2d_symm memory  
  recursive subroutine  SpAMM_delete_tree_2d_symm_recur (self)
    !
    type(tree_2d_symm), pointer, intent(inout) :: self
    type(tree_2d_symm), pointer                :: ch ! sub-tree pointer 
 
    ! check for self-non-association (eg. at leaf pntr) ...
    if(.not.associated(self))return
    
    ch=>self%child_00                          ! take the [00] channel
    call SpAMM_delete_tree_2d_symm_recur (ch)  ! recur to the botom:    
    call SpAMM_delete_tree_2d_node       (ch)  ! kill backwards up the tree 

    ch=>self%child_01                          ! take the [01] channel
    call SpAMM_delete_tree_2d_symm_recur (ch)  ! recur to the botom:    
    call SpAMM_delete_tree_2d_node       (ch)  ! kill backwards up the tree 

    ch=>self%child_11                          ! take the [11] channel
    call SpAMM_delete_tree_2d_symm_recur (ch)  ! recur to the botom:    
    call SpAMM_delete_tree_2d_node       (ch)  ! kill backwards up the tree 

    ! step back up the tree 
  end subroutine SpAMM_delete_tree_2d_symm_recur ! ... and we're out ... 

  subroutine  SpAMM_delete_tree_2d (self)
    type(tree_2d_symm), pointer, intent(inout) :: self
    ! kill the adornment, mort le accoutrement ...    
    if(allocated(self%chunk))deallocate(self%chunk) 
    deallocate(self%frill%bndbx)   ! fru-fru 
    deallocate(self)               ! done 
    nullify(self)                  ! bye-bye
  end subroutine SpAMM_delete_tree_2d


