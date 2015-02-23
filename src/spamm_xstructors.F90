module spamm_xstructors

  use spamm_structures
  use spamm_decoration

  implicit none

contains

  ! SpAMM xstructors ________________ SXTRS _________________
  ! constructors and destructors 
  !
  ! instantiates the tree_2d_symm structure 

  function SpAMM_new_top_tree_2d_symm (NDimn) result (tree)
    !
    integer, dimension(1:2),intent(in) :: NDimn
    integer                            :: M_pad, N_pad, depth
    type(SpAMM_tree_2d_symm), pointer  :: tree

    ! instantiate the root node.  this is the tree top ... 
    allocate(tree)

    ! here are padded dimensions ...
    do depth=0,64
       M_pad=SPAMM_BLOCK_SIZE*2**depth
       if(M_pad>NDimn(1))exit
    enddo
    !
    do depth=0,64
       N_pad=SPAMM_BLOCK_SIZE*2**depth
       if(N_pad>NDimn(2))exit
    enddo

    ! the [i-j] native dimensions ...
    tree%frill%ndimn=ndimn

    ! not a leaf node, this is the top (root) of the tree, k?
    tree%frill%Leaf=.FALSE.

    ! the 2-ary tiles
    tree%frill%bndbx(0:1,1) = (/ 1, M_pad /)  ! [i-lo,i-hi]
    tree%frill%bndbx(0:1,2) = (/ 1, N_pad /)  ! [j-lo,j-hi]

  end function SpAMM_new_top_tree_2d_symm
  
  ! structor for the lo-lo [00] channel ...
  function SpAMM_construct_tree_2d_symm_00(tree) result(ch00)  

    type(SpAMM_tree_2d_symm), target, intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer               :: ch00   
    integer, dimension(:,:),  pointer               :: bb, bb00
    integer, dimension(1:2)                         :: mid,wid

    ch00=>tree%child_00

    if(associated(ch00))return            ! pre-existing?  ok, so later ...
    allocate(ch00)                        ! ... otherwise, instantiate

    ch00%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions

    bb=>tree%frill%bndbx                  ! ... local boxes ... 
    bb00=>ch00%frill%bndbx

    wid(:)=bb(1,:)-bb(0,:)+1              ! ... the 00 split ...
    mid(:)=bb(0,:)+wid(:)/2 

    bb00(:,1)=(/ bb(0,1) , mid(1) /)      ! [lo,mid]
    bb00(:,2)=(/ bb(0,2) , mid(2) /)      ! [lo,mid]
    
    ch00%frill%Leaf=.FALSE.               ! default, not a leaf ...

    ! leaf criterion ... 
    if(mid(1)==SBS)then
       ch00%frill%Leaf=.TRUE.
       allocate(ch00%chunk(1:SBS,1:SBS)) ! grab a chunk for each leaf node, always
       ch00%chunk=SpAMM_Zero        
       ch00%frill%flops=0
    endif

  end function SpAMM_construct_tree_2d_symm_00

  ! structor for the lo-hi [01] channel ...
  function SpAMM_construct_tree_2d_symm_01(tree) result(ch01)  

    type(SpAMM_tree_2d_symm), target, intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer               :: ch01
    integer, dimension(:,:) , pointer               :: bb, bb01
    integer, dimension(1:2)                         :: mid,wid

    ch01=>tree%child_01
    
    if(associated(ch01))return            ! pre-existing?  ok, so later ...
    allocate(ch01)                        ! ... otherwise, instantiate 
    
    ch01%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions

    bb=>tree%frill%bndbx                  ! ... local boxes ... 
    bb01=>ch01%frill%bndbx

    wid(:)=bb(1,:)-bb(0,:)+1              ! ... the 01 split ...
    mid(:)=bb(0,:)+wid(:)/2 

    bb01(:,1)=(/bb(0,1) , mid(1) /)       ! [lo   , mid]
    bb01(:,2)=(/mid(2)+1, bb(1,2)/)       ! [mid+1,  hi]

    ch01%frill%Leaf=.FALSE.               ! default, not a leaf ...

    ! leaf criterion ... 
    if(mid(1)==SBS)then
       ch01%frill%Leaf=.TRUE.
       allocate(ch01%chunk(1:SBS,1:SBS)) ! grab a chunk for the leaf node, always
       ch01%chunk=SpAMM_Zero        
       ch01%frill%flops=0
    endif

  end function SpAMM_construct_tree_2d_symm_01

  ! structor for the hi-hi [11] channel ...
  function SpAMM_construct_tree_2d_symm_11(tree) result(ch11)  

    type(SpAMM_tree_2d_symm), target, intent(inout) :: tree
    type(SpAMM_tree_2d_symm), pointer               :: ch11
    integer, dimension(:,:),  pointer               :: bb, bb11
    integer, dimension(1:2)                         :: mid,wid

    ch11=>tree%child_11
    
    if(associated(ch11))return           ! pre-existing?  ok, so later ...
    allocate(ch11)                       ! ... otherwise, instantiate
    
    ch11%frill%ndimn = tree%frill%ndimn  ! pass down unpadded dimensions

    bb=>tree%frill%bndbx                 ! ... local boxes ... 
    bb11=>ch11%frill%bndbx

    wid(:)=bb(1,:)-bb(0,:)+1             ! ... the 11 split ...
    mid(:)=bb(0,:)+wid(:)/2 

    bb11(0,:)=mid(:)+1                   ! [mid+1, hi]
    bb11(1,:)=bb(1,:)                    ! [mid+1, hi]
    
    ch11%frill%Leaf=.FALSE.              ! default, not a leaf ...

    ! leaf criterion ... 
    if(mid(1)==SBS)then
       ch11%frill%Leaf=.TRUE.
       allocate(ch11%chunk(1:SBS,1:SBS)) ! grab a chunk for the leaf node, always
       ch11%chunk=SpAMM_Zero        
       ch11%frill%flops=0
    endif

  end function SpAMM_construct_tree_2d_symm_11

  ! Structor to recursively destroy a tree_2d_symm ... 
  recursive subroutine  SpAMM_destruct_tree_2d_symm_recur (self)
    !
    type(SpAMM_tree_2d_symm), pointer,  intent(inout) :: self
 
    ! check for self-non-association (eg. at leaf pntr) ...
    if(.not.associated(self))return
    
    call SpAMM_destruct_tree_2d_symm_recur (self%child_00) ! take the [00] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_00) ! kill backwards up the tree

    call SpAMM_destruct_tree_2d_symm_recur (self%child_01) ! take the [01] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_01) ! kill backwards up the tree

    call SpAMM_destruct_tree_2d_symm_recur (self%child_11) ! take the [11] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_11) ! kill backwards up the tree

  end subroutine SpAMM_destruct_tree_2d_symm_recur ! ... and we're out ... 

  ! kill the adornment 
  subroutine  SpAMM_destruct_tree_2d_symm_node (self)

    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: self

    if(self%frill%leaf)deallocate(self%chunk) 
    deallocate(self)               ! done 
    nullify(self)                  ! bye-bye

  end subroutine SpAMM_destruct_tree_2d_symm_node

  ! c => a, b = c => a
  function SpAMM_tree_2d_symm_copy_tree_2d_symm (a, b) result(d)

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)              :: a
    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(INOUT), OPTIONAL :: b
    TYPE(SpAMM_tree_2d_symm), POINTER                          :: d

    IF(PRESENT(b))THEN
       d => b
    ELSEIF(.NOT.ASSOCIATED(a))THEN
       d => null()
       RETURN
    ELSE
       ! nothing passed in, and we have an associated A, so lets pop a new tree top ...
       d => SpAMM_new_top_tree_2d_symm ( b%frill%NDimn )
    ENDIF

    ! d |cpy> a
    CALL SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (d, a)

  END function SpAMM_tree_2d_symm_copy_tree_2d_symm

  !  d |cpy> a
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (d, a)

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: a
    TYPE(SpAMM_tree_2d_symm), POINTER                :: d
    
    if(.not.associated(a))return
    
    if (a%frill%leaf) then       

       d%chunk(1:SBS,1:SBS)=a%chunk(1:SBS,1:SBS) ! d%chunk |cpy> a%chunk

    else

       CALL SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_00(d), a%child_00)
       CALL SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_01(d), a%child_01)
       CALL SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_11(d), a%child_11)

    endif    

    CALL SpAMM_decoration_2d_copy_decoration_2d(d%frill,a%frill) ! d%frill |cpy> a%frill
    
  END SUBROUTINE SpAMM_tree_2d_symm_copy_tree_2d_symm_recur


!!$
!!$  function SpAMM_new_identity_tree_2d_symm ( M, N, Alpha_O ) result (tree)
!!$    !
!!$    integer,        intent(in)  :: M, N
!!$    REAL(SpAMM_KIND), OPTIONAL  :: Alpha_O
!!$    REAL(SpAMM_KIND)            :: Alpha
!!$    integer                     :: depth
!!$    type(SpAMM_tree_2d_symm), pointer :: tree
!!$
!!$    tree => SpAMM_new_top_tree_2d_symm ( M, N )
!!$
!!$    IF(PRESENT(alpha_O))THEN
!!$       alpha=Alpha_O
!!$    ELSE
!!$       alpha=SpAMM_One
!!$    ENDIF    
!!$
!!$    ! push alpha onto the trace ...
!!$    depth=0
!!$    CALL SpAMM_new_identity_tree_2d_symm_recur (tree, alpha, depth)
!!$
!!$  end function SpAMM_new_identity_tree_2d_symm

!!$  ! putting alpha down, onto the trace of this tree_2d_symm ...
!!$  recursive subroutine SpAMM_new_identity_tree_2d_symm_recur (tree, alpha, depth)
!!$
!!$    type(SpAMM_tree_2d_symm)                :: tree
!!$    real(SpAMM_KIND),        intent(in)     :: alpha
!!$    integer,                 intent(in)     :: depth 
!!$    INTEGER, DIMENSION(:,:), pointer        :: bb
!!$    integer                                 :: i 
!!$
!!$    bb => tree%frill%bndbx 
!!$
!!$    IF( bb(1,1)-bb(0,1) == SBS )THEN  ! Leaf condition ? 
!!$
!!$       ! here is a 2d chunk (ch) ... 
!!$       allocate( tree%chunk( SBS, SBS ) ) 
!!$
!!$       ! set its trace with scalar alpha ...
!!$       tree%chunk=SpAMM_Zero
!!$       do i=1, SBS
!!$          tree%chunk(i,i)=alpha 
!!$       enddo
!!$
!!$    ELSEIF(depth>16)THEN ! build the identity tree down
!!$
!!$       ! child along [00]: [lo,mid]x[lo,mid] ... 
!!$       CALL SpAMM_new_identity_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_00(tree) , alpha, depth+1 )
!!$       ! nothing off diagonal 
!!$       tree%child_01=>Null() 
!!$       ! child along [11]: [mid+1,hi]x[mid+1,hi] ... 
!!$       CALL SpAMM_new_identity_tree_2d_symm_recur( SpAMM_construct_tree_2d_symm_11(tree) , alpha, depth+1 )
!!$
!!$    ELSE; STOP ' depth 16 exceeded in SpAMM_init_ident_tree_2d_symm_recur';  
!!$    ENDIF
!!$    
!!$    ! merge & regarnish back up the tree ...
!!$    CALL SpAMM_redecorate_tree_2d_symm(tree)
!!$    !
!!$  end subroutine SpAMM_new_identity_tree_2d_symm_recur
!!$

end module spamm_xstructors

