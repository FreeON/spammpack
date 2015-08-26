!> SpAMM memory operations.
!!
!! Usage example:
!!
!!     a_1 => init (vector top)
module spamm_xstructors

  use spamm_structures
  use spamm_decoration

  implicit none

  !> Occlusion operations.
  interface spamm_occlude
     module procedure spamm_occlude_tree_1d
     module procedure spamm_occlude_tree_2d_symm
     module procedure spamm_occlude_tree_2d_symm_dot_tree_1d
     module procedure spamm_occlude_tree_2d_symm_dot_tree_2d_symm
  end interface spamm_occlude

  !> Flip a tree.
  interface spamm_flip
     module procedure spamm_flip_init_tree_1d_recur
     module procedure spamm_flip_init_tree_2d_symm_recur
  end interface spamm_flip

  !> Prune operations.
  interface spamm_prune
     module procedure spamm_prune_initted_tree_1d_recur
     module procedure spamm_prune_initted_tree_2d_symm_recur
  end interface spamm_prune

contains

  logical function SpAMM_occlude_tree_1d(a, Tau2)

    type(SpAMM_tree_1d), pointer :: a
    real(SPAMM_KIND), intent(IN) :: Tau2

    SpAMM_occlude_tree_1d = .false.
    if(.not. associated(a) )return
    ! cull
    if(a%frill%Norm2 <= Tau2 )return
    ! passed all the checks ...
    SpAMM_occlude_tree_1d = .true.

  end function SpAMM_occlude_tree_1d

  logical function SpAMM_occlude_tree_2d_symm( a, Tau2 )

    type(SpAMM_tree_2d_symm), pointer :: a
    real(SPAMM_KIND),      intent(IN) :: Tau2

    SpAMM_occlude_tree_2d_symm = .false.
    if(.not. associated(a) )return
    ! cull
    if(a%frill%Norm2 <= Tau2 )return
    ! passed all the checks ...
    SpAMM_occlude_tree_2d_symm = .true.

  end function SpAMM_occlude_tree_2d_symm

  logical function SpAMM_occlude_tree_2d_symm_dot_tree_1d( a, b, Tau2 )

    type(SpAMM_tree_2d_symm), pointer, intent(IN) :: a
    type(SpAMM_tree_1d)     , pointer, intent(IN) :: b
    real(SPAMM_KIND),                  intent(IN) :: Tau2

    SpAMM_occlude_tree_2d_symm_dot_tree_1d = .false.

    if(.not. associated(a) )return
    if(.not. associated(b) )return

    ! cull
    if(a%frill%Norm2 * b%frill%Norm2 <= Tau2 )return

    ! passed all the checks ...
    SpAMM_occlude_tree_2d_symm_dot_tree_1d = .true.

  end function SpAMM_occlude_tree_2d_symm_dot_tree_1d

  !> Occlude.
  logical function SpAMM_occlude_tree_2d_symm_dot_tree_2d_symm( a, b, Tau2 )

    type(SpAMM_tree_2d_symm), pointer, intent(IN) :: a,b
    real(SPAMM_KIND),      intent(IN) :: Tau2

    SpAMM_occlude_tree_2d_symm_dot_tree_2d_symm = .false.

    !    write(*,*)associated(a),associated(b)

    if(.not. associated(a)) return
    if(.not. associated(b)) return

    !    write(*,*)a%frill%norm2,b%frill%norm2

    ! cull
    if(a%frill%Norm2*b%frill%Norm2 <= Tau2) return

    ! passed all the checks ...
    SpAMM_occlude_tree_2d_symm_dot_tree_2d_symm = .true.

  end function SpAMM_occlude_tree_2d_symm_dot_tree_2d_symm

  !> Flip on the initialization status of a tree.
  !!
  !! @param a The tree to flip.
  recursive subroutine spamm_flip_init_tree_1d_recur(a)

    type(spamm_tree_1d), pointer :: a

    if(.not.associated(A)) return
    a%frill%needs_initialization = .true.
    call spamm_flip_init_tree_1d_recur(a%child_0)
    call spamm_flip_init_tree_1d_recur(a%child_1)

  end subroutine spamm_flip_init_tree_1d_recur

  !> Flip on the initialization status of a tree.
  !!
  !! @param a The tree to flip.
  recursive subroutine spamm_flip_init_tree_2d_symm_recur(a)

    type(spamm_tree_2d_symm), pointer :: a

    if(.not.associated(a)) return
    a%frill%needs_initialization = .true.
    call spamm_flip_init_tree_2d_symm_recur(a%child_00)
    call spamm_flip_init_tree_2d_symm_recur(a%child_11)
    call spamm_flip_init_tree_2d_symm_recur(a%child_01)
    call spamm_flip_init_tree_2d_symm_recur(a%child_10)

  end subroutine spamm_flip_init_tree_2d_symm_recur

  recursive subroutine SpAMM_Prune_Initted_tree_1d_recur(a)

    type(SpAMM_tree_1d), pointer  :: a

    if(.not.associated(a))return

    if(a%frill%needs_initialization)then
       call SpAMM_destruct_tree_1d_recur (a)
    else
       call SpAMM_Prune_Initted_tree_1d_recur(a%child_0)
       call SpAMM_Prune_Initted_tree_1d_recur(a%child_1)
    end if

  end subroutine SpAMM_Prune_Initted_tree_1d_recur

  recursive subroutine SpAMM_Prune_Initted_tree_2d_symm_recur(a)

    type(SpAMM_tree_2d_symm), pointer  :: a

    if(.not.associated(a))return

    if(a%frill%needs_initialization)then
       call SpAMM_destruct_tree_2d_symm_recur(a)
    else
       call SpAMM_Prune_Initted_tree_2d_symm_recur(a%child_00)
       call SpAMM_Prune_Initted_tree_2d_symm_recur(a%child_11)
       call SpAMM_Prune_Initted_tree_2d_symm_recur(a%child_01)
       call SpAMM_Prune_Initted_tree_2d_symm_recur(a%child_10)
    end if

  end subroutine SpAMM_Prune_Initted_tree_2d_symm_recur

  !> Allocate a new tree.
  !!
  !! @param n The number of elements.
  !! @return The new tree.
  function spamm_new_top_tree_1d(n) result(tree)

    integer, intent(in) :: n
    type(spamm_tree_1d), pointer :: tree

    integer :: n_pad, depth

    ! Instantiate the root node.  This is the tree top ...
    allocate(tree)

    ! Here are padded dimensions ...
    do depth = 0, 64
       n_pad = SPAMM_CHUNK_SIZE*2**depth
       if(n_pad >= n) exit
    end do

    ! The [i] native dimension ...
    tree%frill%ndimn = n

    ! The [i] padded width
    tree%frill%width = n_pad

    ! The native tile
    tree%frill%bndbx(0:1) = [1, n]  ! [i-lo,i-hi]

  end function spamm_new_top_tree_1d

  !++XSTRUCTORS:     SpAMM_construct_tree_1d_0
  !++XSTRUCTORS:       a_1%0 => init (channel [0] constructor)
  function spamm_construct_child_1d(tree, child_index) result(child)

    type(spamm_tree_1d), pointer :: tree
    integer, intent(in) :: child_index
    type(spamm_tree_1d), pointer :: child

    integer :: mi

    select case(child_index)
    case(0)
       if(associated(tree%child_0)) then
          child => tree%child_0
          return ! pre-existing?  ok, so later ...
       end if
    case(1)
       if(associated(tree%child_1)) then
          child => tree%child_1
          return ! pre-existing?  ok, so later ...
       end if
    end select

    associate(lo => tree%frill%bndbx(0), &
         hi => tree%frill%bndbx(1), &
         wi => tree%frill%width)

      mi = min(lo+wi/2-1, hi)
      if(mi+1 > hi) then
         child => null()
         return                                    ! margin over-run
      end if

      allocate(child)

      child%frill%width = wi/2
      child%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions
      select case(child_index)
      case(0)
         child%frill%bndbx = [lo, mi]
      case(1)
         child%frill%bndbx = [mi+1, hi]
      end select

      child%frill%leaf = .false.
      child%frill%flops = SPAMM_INIT
      child%frill%norm2 = SPAMM_INIT
      if(wi == 2*SBS) then                          ! leaf criterion ...
         child%frill%leaf = .true.
         allocate(child%chunk(SBS))          ! grab a chunk for each leaf node, always
         child%chunk = 0
      end if

      select case(child_index)
      case(0)
         tree%child_0 => child
      case(1)
         tree%child_1 => child
      end select
    end associate

  end function spamm_construct_child_1d

  recursive subroutine SpAMM_destruct_tree_1d_recur(self)   !++
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

    if(.not.associated(self))return
    if(self%frill%leaf)deallocate(self%chunk)
    deallocate(self)
    nullify(self)                  ! bye-bye

  end subroutine SpAMM_destruct_tree_1d_node

  !++XSTRUCTORS:     SpAMM_tree_1d_copy_tree_1d
  !++XSTRUCTORS:       d_1 => a (wrapper)
  function SpAMM_tree_1d_copy_tree_1d (a, in_O) result(d)

    type(SpAMM_tree_1d), pointer, intent(IN)              :: a
    type(SpAMM_tree_1d), pointer, intent(INOUT), optional :: in_O
    type(SpAMM_tree_1d), pointer                          :: d

    d => null()
    if(present(in_O))then
       d => in_O
    elseif(.not.associated(a))then
       return
    end if

    ! nothing passed in, and we have an associated A, so lets pop a new tree top ...
    if(.not.associated(d)) &
         d => SpAMM_new_top_tree_1d ( a%frill%NDimn )

    ! d |cpy> a
    call SpAMM_tree_1d_copy_tree_1d_recur (d, a)

  end function SpAMM_tree_1d_copy_tree_1d

  !++XSTRUCTORS:     SpAMM_tree_1d_copy_tree_1d_recur
  !++XSTRUCTORS:       d_1 => a (recursive)
  recursive subroutine SpAMM_tree_1d_copy_tree_1d_recur (d, a)

    type(SpAMM_tree_1d), pointer, intent(IN)    :: a
    type(SpAMM_tree_1d), pointer                :: d

    if(.not.associated(a))return

    if(a%frill%leaf ) then

       d%chunk(1:SBS)=a%chunk(1:SBS)

    else

       !       IF(ASSOCIATED(a%child_0))&
       call spamm_tree_1d_copy_tree_1d_recur(spamm_construct_child_1d(d, 0), a%child_0)
       !       IF(ASSOCIATED(a%child_1))&
       call spamm_tree_1d_copy_tree_1d_recur(spamm_construct_child_1d(d, 1), a%child_1)

    end if

    call SpAMM_redecorate_tree_1d(d)

  end subroutine SpAMM_tree_1d_copy_tree_1d_recur

  !> Construct new 2D tree.
  !!
  !! @param ndimn The matrix size.
  !! @return The new matrix.
  function spamm_new_top_tree_2d_symm(ndimn) result(tree)

    integer, intent(in) :: ndimn(2)
    type(spamm_tree_2d_symm), pointer :: tree

    integer :: m_pad, n_pad, depth

    ! instantiate the root node. This is the tree top ...
    allocate(tree)

    ! here are padded dimensions ...
    do depth = 0, 64
       m_pad = SPAMM_CHUNK_SIZE*2**depth
       if(m_pad >= ndimn(1)) exit
    end do

    do depth = 0, 64
       n_pad = SPAMM_CHUNK_SIZE*2**depth
       if(n_pad >= ndimn(2)) exit
    end do

    ! the [i]-[j] native dimensions ...
    tree%frill%ndimn = ndimn

    ! the [i]-[j] padded width
    tree%frill%width = [M_pad, N_pad]

    ! the 2-ary tiles
    tree%frill%bndbx(0, :) = [1, 1]  ! [lo,lo]
    tree%frill%bndbx(1, :) = NDimn   ! [hi,hi]

    ! inited measures
    tree%frill%non0s = SpAMM_init
    tree%frill%norm2 = SpAMM_init
    tree%frill%flops = SpAMM_init

    ! check that we might the top may be the leaf
    if(SBS >= NDimn(1)) then
       tree%frill%leaf = .true.
       allocate(tree%chunk(1:SBS, 1:SBS))  ! leaf == allocated(chunk)
       tree%chunk = 0              ! init
    end if

    ! no kids
    tree%child_00 => null()
    tree%child_01 => null()
    tree%child_10 => null()
    tree%child_11 => null()

  end function spamm_new_top_tree_2d_symm

  !> SpAMM_construct_tree_2d_symm_00
  !! a_2%00 => init (constructor of the lo-lo [00] channel)
  function spamm_construct_tree_2d_symm_00(tree) result(ch00)

    type(spamm_tree_2d_symm), pointer :: tree
    type(spamm_tree_2d_symm), pointer :: ch00

    integer, dimension(1:2) :: lo, hi, mi, wi
    integer :: i

    if(associated(tree%child_00)) then
       ch00 => tree%child_00
       return                                      ! pre-existing?  ok, so later ...
    end if
    allocate(tree%child_00)                        ! ... otherwise, instantiate
    lo = tree%frill%bndbx(0, :)
    hi = tree%frill%bndbx(1, :)
    wi = tree%frill%width
    mi = lo+wi/2-1
    do i = 1, 2
       mi(i) = min(hi(i), mi(i))
    end do
    tree%child_00%frill%needs_initialization = .true. ! a new node, so set init status true ...
    tree%child_00%frill%width = wi/2               ! next level width
    tree%child_00%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions
    tree%child_00%frill%bndbx(:, 1) = (/lo(1), mi(1)/) ! [lo:mid][i]
    tree%child_00%frill%bndbx(:, 2) = (/lo(2), mi(2)/) ! [lo:mid][j]
    tree%child_00%frill%leaf = .false.               ! default, not a leaf
    tree%child_00%frill%flops = SPAMM_INIT
    tree%child_00%frill%norm2 = SPAMM_INIT
    if(wi(1) == 2*SBS) then                           ! at resolution?
       tree%child_00%frill%leaf = .true.             ! we have a leaf
       allocate(tree%child_00%chunk(1:SBS, 1:SBS))  ! leaf == allocated(chunk)
       tree%child_00%chunk = 0              ! init
    end if
    ch00 => tree%child_00

  end function spamm_construct_tree_2d_symm_00

  !++XSTRUCTORS:     SpAMM_construct_tree_2d_symm_01
  !++XSTRUCTORS:       a_2%01 => init (constructor of the lo-hi [01] channel)
  function spamm_construct_tree_2d_symm_01(tree) result(ch01)

    type(SpAMM_tree_2d_symm), pointer :: tree
    type(SpAMM_tree_2d_symm), pointer :: ch01
    integer, dimension(1:2)           :: lo,hi,mi,wi

    if(associated(tree%child_01))then
       ch01=>tree%child_01
       return                                      ! pre-existing?  ok, so later ...
    end if

    lo = tree%frill%bndbx(0,:)
    hi = tree%frill%bndbx(1,:)

    !    if(hi(1)>876.or.hi(2)>876)then
    !       stop '01'
    !   end if

    wi = tree%frill%width
    mi = lo+wi/2-1

    mi(1)=min(hi(1),mi(1))
    if(mi(2)+1>hi(2))then
       ch01=>NULL()
       return                                     ! margin over-run
    end if
    allocate(tree%child_01)                       ! ... otherwise, instantiate

    tree%child_01%frill%needs_initialization = .true.              ! a new node, so set init status true ...
    tree%child_01%frill%width = wi/2               ! next level width
    tree%child_01%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions
    tree%child_01%frill%bndbx(:,1)=(/lo(1)  ,mi(1)/) ! [lo   ,mid][i]
    tree%child_01%frill%bndbx(:,2)=(/mi(2)+1,hi(2)/) ! [mid+1, hi][j]
    tree%child_01%frill%Leaf=.false.               ! default, not a leaf
    tree%child_01%frill%flops=SpAMM_init
    tree%child_01%frill%norm2=SpAMM_init
    if(wi(1)==2*SBS)then                           ! at resolution?
       tree%child_01%frill%Leaf=.true.             ! we have a leaf
       allocate(tree%child_01%chunk(SBS,SBS))      ! leaf == allocated(chunk)
       tree%child_01%chunk=0                       ! init
    end if

    !    write(*,33) tree%child_01%frill%bndbx(:,1),tree%child_01%frill%bndbx(:,2), &
    !            wi,tree%child_01%frill%leaf
    !33  format(' 01: [ ',I3,", ",I3," ]x[ ",I3,", ",I3," ], wid = ",2I4,4L3 )

    ch01=>tree%child_01

  end function SpAMM_construct_tree_2d_symm_01

  function SpAMM_construct_tree_2d_symm_10(tree) result(ch10)

    type(SpAMM_tree_2d_symm), pointer :: tree
    type(SpAMM_tree_2d_symm), pointer :: ch10
    integer, dimension(1:2)           :: lo,hi,mi,wi

    if(associated(tree%child_10))then
       ch10=>tree%child_10
       return                                      ! pre-existing?  ok, so later ...
    end if

    lo = tree%frill%bndbx(0,:)
    hi = tree%frill%bndbx(1,:)

    wi = tree%frill%width
    mi = lo+wi/2-1

    mi(2)=min(hi(2),mi(2))
    if(mi(1)+1>hi(1))then
       ch10=>NULL()
       return                                     ! margin over-run
    end if
    allocate(tree%child_10)                       ! ... otherwise, instantiate

    tree%child_10%frill%needs_initialization = .true.              ! a new node, so set init status true ...
    tree%child_10%frill%width = wi/2               ! next level width
    tree%child_10%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions
    tree%child_10%frill%bndbx(:,1)=(/mi(1)+1,hi(1)/) ! [mid+1, hi][i]
    tree%child_10%frill%bndbx(:,2)=(/lo(2)  ,mi(2)/) ! [lo   ,mid][j]
    tree%child_10%frill%Leaf=.false.               ! default, not a leaf
    tree%child_10%frill%flops=SpAMM_init
    tree%child_10%frill%norm2=SpAMM_init
    if(wi(1)==2*SBS)then                           ! at resolution?
       tree%child_10%frill%Leaf=.true.             ! we have a leaf
       allocate(tree%child_10%chunk(1:SBS,1:SBS))  ! leaf == allocated(chunk)
       tree%child_10%chunk=0              ! init
    end if

    !    write(*,33) tree%child_10%frill%bndbx(:,1),tree%child_10%frill%bndbx(:,2), &
    !            wi,tree%child_10%frill%leaf
    !33  format(' 10: [ ',I3,", ",I3," ]x[ ",I3,", ",I3," ], wid = ",2I4,4L3 )

    ch10=>tree%child_10

  end function SpAMM_construct_tree_2d_symm_10

  !++XSTRUCTORS:     SpAMM_construct_tree_2d_symm_01
  !++XSTRUCTORS:       a_2%11 => init (constructor of the hi-hi [11] channel)
  function SpAMM_construct_tree_2d_symm_11(tree) result(ch11)

    type(SpAMM_tree_2d_symm), pointer :: tree
    type(SpAMM_tree_2d_symm), pointer :: ch11
    integer, dimension(1:2)           :: lo,hi,mi,wi

    if(associated(tree%child_11))then
       ch11=>tree%child_11
       return                                     ! pre-existing?  ok, so later ...
    end if

    lo = tree%frill%bndbx(0,:)
    hi = tree%frill%bndbx(1,:)

    wi = tree%frill%width
    mi = lo+wi/2-1

    mi(1)=min(hi(1),mi(1))
    mi(2)=min(hi(2),mi(2))
    if(mi(1)+1>hi(1) .or. mi(2)+1>hi(2))then
       ch11=>NULL()
       return                                    ! margin over-run
    end if
    allocate(tree%child_11)                       ! ... otherwise, instantiate

    tree%child_11%frill%needs_initialization = .true.              ! a new node, so set init status true ...
    tree%child_11%frill%width = wi/2               ! next level width
    tree%child_11%frill%ndimn = tree%frill%ndimn   ! pass down unpadded dimensions
    tree%child_11%frill%bndbx(0,:)=mi(:)+1                ! [mid+1, hi]
    tree%child_11%frill%bndbx(1,:)=hi                      ! [mid+1, hi]
    tree%child_11%frill%Leaf=.false.               ! default, not a leaf
    tree%child_11%frill%flops=SpAMM_init
    tree%child_11%frill%norm2=SpAMM_init
    if(wi(1)==2*SBS)then                           ! at resolution?
       tree%child_11%frill%Leaf=.true.             ! we have a leaf
       allocate(tree%child_11%chunk(1:SBS,1:SBS))  ! leaf == allocated(chunk)
       tree%child_11%chunk=0              ! init
    end if

    !    write(*,33) tree%child_11%frill%bndbx(:,1) ,tree%child_11%frill%bndbx(:,2),wi/2,tree%child_11%frill%leaf
    !33  format(' 11: [ ',I3,", ",I3," ]x[ ",I3,", ",I3," ], wid = ",2I4,4L3 )

    ch11=>tree%child_11

  end function SpAMM_construct_tree_2d_symm_11

  !++XSTRUCTORS:     SpAMM_destruct_tree_2d_symm_recur
  !++XSTRUCTORS:       a_2 => null() (recursive destructor of the symmetric matrix)
  recursive subroutine spamm_destruct_tree_2d_symm_recur(self)

    type(spamm_tree_2d_symm), pointer, intent(inout) :: self

    ! check for self-non-association (eg. at leaf pntr) ...
    if(.not. associated(self)) return

    call SpAMM_destruct_tree_2d_symm_recur (self%child_00) ! take the [00] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_00) ! kill backwards up the tree

    call SpAMM_destruct_tree_2d_symm_recur (self%child_01) ! take the [01] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_01) ! kill backwards up the tree

    call SpAMM_destruct_tree_2d_symm_recur (self%child_10) ! take the [10] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_10) ! kill backwards up the tree

    call SpAMM_destruct_tree_2d_symm_recur (self%child_11) ! take the [11] channel
    call SpAMM_destruct_tree_2d_symm_node  (self%child_11) ! kill backwards up the tree

    call SpAMM_destruct_tree_2d_symm_node  (self)          ! kill y' self

  end subroutine SpAMM_destruct_tree_2d_symm_recur ! ... and we're out ...

  !++XSTRUCTORS:     SpAMM_destruct_tree_2d_symm_node
  !++XSTRUCTORS:       a_2 => null() (node level destructor of the symmetric matrix)
  subroutine  SpAMM_destruct_tree_2d_symm_node (self)

    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: self

    if(.not.associated(self))return
    if(self%frill%leaf)deallocate(self%chunk)
    deallocate(self)               ! done
    nullify(self)                  ! bye-bye

  end subroutine SpAMM_destruct_tree_2d_symm_node

  ! - - - - - - - - - - - - - - - - -2d, 2d, 2d - - - - - - - - - - - - - - - - - -


  !++XSTRUCTORS:     SpAMM_tree_2d_symm_copy_tree_2d_symm
  !++XSTRUCTORS:       d_2 => a_2  (wrapper)
  function SpAMM_tree_2d_symm_copy_tree_2d_symm (a, in_O, threshold_O, symmetrize_O ) result(d)

    type(SpAMM_tree_2d_symm), pointer, intent(IN)              :: a
    real(SPAMM_KIND),                  intent(IN),    optional :: threshold_o
    logical,                           intent(IN),    optional :: symmetrize_O
    type(SpAMM_tree_2d_symm), pointer, intent(INOUT), optional :: in_O
    type(SpAMM_tree_2d_symm), pointer                          :: d
    real(SPAMM_KIND)                                           :: threshold2

    d => null()
    if(present(in_O))then
       d => in_O
    elseif(.not.associated(a))then
       return
    end if

    if(.not.associated(d)) d => SpAMM_new_top_tree_2d_symm (a%frill%NDimn )

    ! d |cpy>|threshold?> a
    threshold2=0
    if(present(threshold_o))threshold2=threshold_O**2

    call SpAMM_flip(d)

    if(present(symmetrize_O))then
       if(symmetrize_O)then
          stop ' not yet'
          call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize (d, a, threshold2)
       else
          call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (d, a, threshold2)
       end if
    else
       call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (d, a, threshold2)
    end if

    call SpAMM_prune(d)

  end function SpAMM_tree_2d_symm_copy_tree_2d_symm

  !++XSTRUCTORS:     SpAMM_tree_2d_symm_copy_tree_2d_symm_recur
  !++XSTRUCTORS:       d_2 => a_2  (recursive)
  recursive subroutine SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (d, a, Tau2)

    type(SpAMM_tree_2d_symm), pointer, intent(IN)    :: a
    real(SPAMM_KIND),                  intent(IN)    :: Tau2
    type(SpAMM_tree_2d_symm), pointer                :: a00,a11,a01,a10
    type(SpAMM_tree_2d_symm), pointer                :: d
    type(SpAMM_tree_2d_symm), pointer                :: d00,d11,d01,d10

    if(a%frill%leaf)then
       d%frill%needs_initialization=.false.
       d%chunk(1:SBS,1:SBS)=a%chunk(1:SBS,1:SBS) ! d%chunk |cpy> a%chunk
       d%frill%flops=0
    else
       ! local children
       a00=>a%child_00; a11=>a%child_11; a01=>a%child_01; a10=>a%child_10;
       d00=>NULL();     d11=>NULL();     d01=>NULL();     d10=>NULL();
       ! copy diagonal
       if(SpAMM_occlude( a00, Tau2 ) ) &
            call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_00(d), a00, Tau2 )
       if(SpAMM_occlude( a11, Tau2 ) ) &
            call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_11(d), a11, Tau2 )
       if(SpAMM_occlude( a01, Tau2 ) ) &
            call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_01(d), a01, Tau2 )
       if(SpAMM_occlude( a10, Tau2 ) ) &
            call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur (SpAMM_construct_tree_2d_symm_10(d), a10, Tau2 )
    end if

    ! redecorate
    call SpAMM_redecorate_tree_2d_symm(d)

  end subroutine SpAMM_tree_2d_symm_copy_tree_2d_symm_recur


  !++XSTRUCTORS:     SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize
  !++XSTRUCTORS:       d_2 => a_2  (recursive)
  recursive subroutine SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize (d, a, threshold2, at)

    type(SpAMM_tree_2d_symm), pointer, intent(IN)           :: a
    type(SpAMM_tree_2d_symm), pointer, intent(IN), optional :: at
    real(SPAMM_KIND),                  intent(IN)    :: threshold2
    type(SpAMM_tree_2d_symm), pointer                :: d

    if(.not.associated(a).and..not.associated(d))then
       return
    elseif(.not.associated(a).and.associated(d))then
       call SpAMM_destruct_tree_2d_symm_recur (d)
       return
    elseif(a%frill%norm2 <= threshold2)then
       call SpAMM_destruct_tree_2d_symm_recur (d)
       return
    elseif(a%frill%leaf) then
       d%chunk(1:SBS,1:SBS)=(a%chunk(1:SBS,1:SBS)+transpose(  a%chunk(1:SBS,1:SBS) ))*SpAMM_half
       ! flops
    else
       call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize (SpAMM_construct_tree_2d_symm_00(d), &
            a%child_00, threshold2 )
       call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize (SpAMM_construct_tree_2d_symm_10(d), &
            a%child_10, threshold2, at=a%child_01 )
       call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize (SpAMM_construct_tree_2d_symm_01(d), &
            a%child_01, threshold2, at=a%child_10 )
       call SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize (SpAMM_construct_tree_2d_symm_11(d), &
            a%child_11, threshold2 )
    end if

    call SpAMM_redecorate_tree_2d_symm(d)

  end subroutine SpAMM_tree_2d_symm_copy_tree_2d_symm_recur_symmetrize


end module spamm_xstructors
