module spamm_xstructors

  use spamm_structures
  use spamm_decoration

  implicit none

contains
  !++XSTRUCTORS: SpAMM memory opperations _______________________ XSTRUCTORS _________________
  !++XSTRUCTORS: constructors and destructors for SpAMM tree-nd structures ...
  !++XSTRUCTORS:   ... TREE-ONE-D ... TREE-ONE-D ... TREE-ONE-D ... 
  !++XSTRUCTORS:     SpAMM_new_top_tree_1d 
  !++XSTRUCTORS:       a_1 => init (vector top) 
  function SpAMM_new_top_tree_1d(NDimn) result (tree)
    !
    integer                       :: NDimn
    integer                       :: M_pad, depth
    type(SpAMM_tree_1d), pointer  :: tree

    ! instantiate the root node.  this is the tree top ... 
    allocate(tree)

    ! here are padded dimensions ...
    do depth=0,64
       M_pad=SPAMM_BLOCK_SIZE*2**depth
       if(M_pad>NDimn)exit
    enddo

    ! the [i] native dimension ...
    tree%frill%ndimn=ndimn

    ! not a leaf node, this is the top (root) of the tree, k?
    tree%frill%Leaf=.FALSE.

    ! the 1-ary tile
    tree%frill%bndbx(0:1)=(/ 1, M_pad /)  ! [i-lo,i-hi]

    ! no kids yet
    tree%child_0=>NULL()
    tree%child_1=>NULL()

  end function SpAMM_new_top_tree_1d

  !++XSTRUCTORS:     SpAMM_construct_tree_1d_0
  !++XSTRUCTORS:       a_1%0 => init (constructor of the lo [0] channel)
  function SpAMM_construct_tree_1d_0(tree) result(ch0)  

    type(SpAMM_tree_1d),    pointer :: tree
    type(SpAMM_tree_1d),    pointer :: ch0   
    integer                         :: mid,wid


    if(associated(tree%child_0))then
       ch0=>tree%child_0
       return            ! pre-existing?  ok, so later ...
    endif
    allocate(tree%child_0)                        ! ... otherwise, instantiate

    tree%child_0%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions

    wid=tree%frill%bndbx(1)-tree%frill%bndbx(0)+1                    ! ... the 0 split ...
    mid=tree%frill%bndbx(0)+wid/2-1 

    tree%child_0%frill%bndbx(:)=(/tree%frill%bndbx(0),mid/)                 ! [lo,mid]
    
    tree%child_0%frill%Leaf=.FALSE.               ! default, not a leaf ...

    ! leaf criterion ... 
    if(wid==2*SBS)then
       tree%child_0%frill%Leaf=.TRUE.
       allocate(tree%child_0%chunk(1:SBS))        ! grab a chunk for each leaf node, always
       tree%child_0%chunk=SpAMM_Zero        
       tree%child_0%frill%flops=0
    endif

    ch0=>tree%child_0

  end function SpAMM_construct_tree_1d_0

  !++XSTRUCTORS:     SpAMM_construct_tree_1d_1
  !++XSTRUCTORS:       a_2%1 => init (constructor of the hi [1] channel)
  function SpAMM_construct_tree_1d_1(tree) result(ch1)  

    type(SpAMM_tree_1d), pointer :: tree
    type(SpAMM_tree_1d), pointer :: ch1
    integer                      :: mid,wid
    
    if(associated(tree%child_1))then       
       ch1=>tree%child_1
       return           ! pre-existing?  ok, so later ...
    endif
    allocate(tree%child_1)                       ! ... otherwise, instantiate
    
    tree%child_1%frill%ndimn = tree%frill%ndimn  ! pass down unpadded dimensions

    wid=tree%frill%bndbx(1)-tree%frill%bndbx(0)+1                   ! ... the 1 split ...
    mid=tree%frill%bndbx(0)+wid/2-1 
    tree%child_1%frill%bndbx(:)=(/mid+1, tree%frill%bndbx(1)/)             ! [mid+1, hi]
    
    tree%child_1%frill%Leaf=.FALSE.              ! default, not a leaf ...

    ! leaf criterion ... 
    if(wid==2*SBS)then
       tree%child_1%frill%Leaf=.TRUE.
       allocate(tree%child_1%chunk(1:SBS))       ! grab a chunk for the leaf node, always
       tree%child_1%chunk=SpAMM_Zero        
       tree%child_1%frill%flops=0
    endif

    ch1=>tree%child_1

  end function SpAMM_construct_tree_1d_1

  recursive subroutine  SpAMM_destruct_tree_1d_recur (self)   !++
  !++XSTRUCTORS:       a_1 => null() (recursive vector destruction )  
    !
    type(SpAMM_tree_1d), pointer,  intent(inout) :: self
 
    ! check for self-non-association (eg. at leaf pntr) ...
    if(.not.associated(self))return
    
    call SpAMM_destruct_tree_1d_recur (self%child_0) ! take the [0] channel
    call SpAMM_destruct_tree_1d_node  (self%child_0) ! kill backwards up the tree

    call SpAMM_destruct_tree_1d_recur (self%child_1) ! take the [1] channel
    call SpAMM_destruct_tree_1d_node  (self%child_1) ! kill backwards up the tree

  end subroutine SpAMM_destruct_tree_1d_recur ! ... and we're out ... 

  subroutine  SpAMM_destruct_tree_1d_node (self)    !++
  !++XSTRUCTORS:       a_1 => null() (node level destructor of the symmetric matrix)  

    type(SpAMM_tree_1d), pointer, intent(inout) :: self

    if(self%frill%leaf)deallocate(self%chunk) 
    deallocate(self)               ! done 
    nullify(self)                  ! bye-bye

  end subroutine SpAMM_destruct_tree_1d_node

  !++XSTRUCTORS:     SpAMM_tree_1d_copy_tree_1d
  !++XSTRUCTORS:       d_1 => a (wrapper)
  function SpAMM_tree_1d_copy_tree_1d (a, b) result(d)

    TYPE(SpAMM_tree_1d), POINTER, INTENT(IN)              :: a
    TYPE(SpAMM_tree_1d), POINTER, INTENT(INOUT), OPTIONAL :: b
    TYPE(SpAMM_tree_1d), POINTER                          :: d

    IF(PRESENT(b))THEN
       d => b
    ELSEIF(.NOT.ASSOCIATED(a))THEN
       d => null()
       RETURN
    ELSE
       ! nothing passed in, and we have an associated A, so lets pop a new tree top ...
       d => SpAMM_new_top_tree_1d ( b%frill%NDimn )
    ENDIF

    ! d |cpy> a
    CALL SpAMM_tree_1d_copy_tree_1d_recur (d, a)

  END function SpAMM_tree_1d_copy_tree_1d

  !++XSTRUCTORS:     SpAMM_tree_1d_copy_tree_1d_recur 
  !++XSTRUCTORS:       d_1 => a (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_1d_copy_tree_1d_recur (d, a)

    TYPE(SpAMM_tree_1d), POINTER, INTENT(IN)    :: a
    TYPE(SpAMM_tree_1d), POINTER                :: d
    
    if(.not.associated(a))return
    
    if (a%frill%leaf) then       
       d%chunk(1:SBS)=a%chunk(1:SBS) 
    else

       CALL SpAMM_tree_1d_copy_tree_1d_recur (SpAMM_construct_tree_1d_0(d), a%child_0)
       CALL SpAMM_tree_1d_copy_tree_1d_recur (SpAMM_construct_tree_1d_1(d), a%child_1)

    endif    

    CALL SpAMM_decoration_1d_copy_decoration_1d(d%frill,a%frill) ! d%frill |cpy> a%frill
    
  END SUBROUTINE SpAMM_tree_1d_copy_tree_1d_recur

  !!
  !++XSTRUCTORS:   ... TREE-TWO-D ... TREE-TWO-D ... TREE-TWO-D ... 
  !++XSTRUCTORS:     SpAMM_new_top_tree_2d_symm 
  !++XSTRUCTORS:       a_2 => init (matrix top)
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

    ! no kids yet
    tree%child_00=>NULL()
    tree%child_01=>NULL()
    tree%child_11=>NULL()

  end function SpAMM_new_top_tree_2d_symm

  !++XSTRUCTORS:     SpAMM_construct_tree_2d_symm_00
  !++XSTRUCTORS:       a_2%00 => init (constructor of the lo-lo [00] channel)
  function SpAMM_construct_tree_2d_symm_00(tree) result(ch00)  

    type(SpAMM_tree_2d_symm), POINTER  :: tree
    type(SpAMM_tree_2d_symm), POINTER  :: ch00
    integer, dimension(1:2)            :: mid,wid

    if(associated(tree%child_00))return            ! pre-existing?  ok, so later ...
    allocate(tree%child_00)                        ! ... otherwise, instantiate

    tree%child_00%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions

    wid(:)=tree%frill%bndbx(1,:)-tree%frill%bndbx(0,:)+1              ! ... the 00 split ...
    mid(:)=tree%frill%bndbx(0,:)+wid(:)/2-1 

    tree%child_00%frill%bndbx(:,1)=(/ tree%frill%bndbx(0,1) , mid(1) /)      ! [lo,mid]
    tree%child_00%frill%bndbx(:,2)=(/ tree%frill%bndbx(0,2) , mid(2) /)      ! [lo,mid]
    
    tree%child_00%frill%Leaf=.FALSE.               ! default, not a leaf ...

    ! leaf criterion ... 
    if(wid(1)==2*SBS)then
       tree%child_00%frill%Leaf=.TRUE.
       allocate(tree%child_00%chunk(1:SBS,1:SBS))  ! grab a chunk for each leaf node, always
       tree%child_00%chunk=SpAMM_Zero        
       tree%child_00%frill%flops=0
    endif

    ch00=>tree%child_00

  end function SpAMM_construct_tree_2d_symm_00
 

  !++XSTRUCTORS:     SpAMM_construct_tree_2d_symm_01
  !++XSTRUCTORS:       a_2%01 => init (constructor of the lo-hi [01] channel)
  function SpAMM_construct_tree_2d_symm_01(tree) result(ch01)  

    type(SpAMM_tree_2d_symm), pointer  :: tree
    type(SpAMM_tree_2d_symm), pointer  :: ch01
    integer, dimension(1:2)            :: mid,wid

    if(associated(tree%child_01))return            ! pre-existing?  ok, so later ...
    allocate(tree%child_01)                        ! ... otherwise, instantiate 
    
    tree%child_01%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions

    wid(:)=tree%frill%bndbx(1,:)-tree%frill%bndbx(0,:)+1  ! ... the 01 split ...
    mid(:)=tree%frill%bndbx(0,:)+wid(:)/2-1 

    tree%child_01%frill%bndbx(:,1)=(/tree%frill%bndbx(0,1) , mid(1) /)       ! [lo   , mid]
    tree%child_01%frill%bndbx(:,2)=(/mid(2)+1, tree%frill%bndbx(1,2)/)       ! [mid+1,  hi]

    tree%child_01%frill%Leaf=.FALSE.               ! default, not a leaf ...

    ! leaf criterion ... 
    if(wid(1)==2*SBS)then
       tree%child_01%frill%Leaf=.TRUE.
       allocate(tree%child_01%chunk(1:SBS,1:SBS))  ! grab a chunk for the leaf node, always
       tree%child_01%chunk=SpAMM_Zero        
       tree%child_01%frill%flops=0
    endif

    ch01=>tree%child_01
    
  end function SpAMM_construct_tree_2d_symm_01

  !++XSTRUCTORS:     SpAMM_construct_tree_2d_symm_01
  !++XSTRUCTORS:       a_2%11 => init (constructor of the hi-hi [11] channel)
  function SpAMM_construct_tree_2d_symm_11(tree) result(ch11)  

    type(SpAMM_tree_2d_symm), pointer :: tree
    type(SpAMM_tree_2d_symm), pointer :: ch11
    integer, dimension(1:2)           :: mid,wid

    if(associated(tree%child_11))return           ! pre-existing?  ok, so later ...
    allocate(tree%child_11)                       ! ... otherwise, instantiate
    
    tree%child_11%frill%ndimn = tree%frill%ndimn  ! pass down unpadded dimensions

    wid(:)=tree%frill%bndbx(1,:)-tree%frill%bndbx(0,:)+1             ! ... the 11 split ...
    mid(:)=tree%frill%bndbx(0,:)+wid(:)/2-1 

    tree%child_11%frill%bndbx(0,:)=mid(:)+1                   ! [mid+1, hi]
    tree%child_11%frill%bndbx(1,:)=tree%frill%bndbx(1,:)                    ! [mid+1, hi]
    
    tree%child_11%frill%Leaf=.FALSE.              ! default, not a leaf ...

    ! leaf criterion ... 
    if(wid(1)==2*SBS)then
       tree%child_11%frill%Leaf=.TRUE.
       allocate(tree%child_11%chunk(1:SBS,1:SBS)) ! grab a chunk for the leaf node, always
       tree%child_11%chunk=SpAMM_Zero        
       tree%child_11%frill%flops=0
    endif

    ch11=>tree%child_11

  end function SpAMM_construct_tree_2d_symm_11

  !++XSTRUCTORS:     SpAMM_destruct_tree_2d_symm_recur 
  !++XSTRUCTORS:       a_2 => null() (recursive destructor of the symmetric matrix)  
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

  !++XSTRUCTORS:     SpAMM_destruct_tree_2d_symm_node
  !++XSTRUCTORS:       a_2 => null() (node level destructor of the symmetric matrix)  
  subroutine  SpAMM_destruct_tree_2d_symm_node (self)

    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: self

    if(self%frill%leaf)deallocate(self%chunk) 
    deallocate(self)               ! done 
    nullify(self)                  ! bye-bye

  end subroutine SpAMM_destruct_tree_2d_symm_node

  !++XSTRUCTORS:     SpAMM_tree_2d_symm_copy_tree_2d_symm
  !++XSTRUCTORS:       d_2 => a_2  (wrapper)
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

  !++XSTRUCTORS:     SpAMM_tree_2d_symm_copy_tree_2d_symm_recur 
  !++XSTRUCTORS:       d_2 => a_2  (recursive)
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

