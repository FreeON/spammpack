#define SpAMM_PRINT_STREAM
module spamm_nbdyalgbra_times

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals

  implicit none


#ifdef SpAMM_PRINT_STREAM
  !> Globals on only in this non-production instance ...
  type SpAMM_cubes
     !> Next element.
     type(SpAMM_cubes), pointer  :: Next
     !> The depth in the tree.
     integer                     :: Levl
     real(SPAMM_KIND)            :: Size
     integer,     dimension(3)   :: Lw, Hi
  end type SpAMM_cubes

  type(SpAMM_cubes), pointer     :: Stream
  integer :: number_stream_elements
  logical :: Build_Stream

#endif

contains

  !++NBODYTIMES: SpAMM generalized n-body algebras for times ____ NBODYTIMES _________________
  !++NBODYTIMES: generalized products, dots, contractions & convolutions (X)
  !++NBODYTIMES:   ... [TREE-ONE-D X TREE-ONE-D] ... [TREE-ONE-D X TREE-ONE-D] ...
  !++NBODYTIMES:   SpAMM_tree_1d_dot_tree_1d_recur
  !++NBODYTIMES:     dot = (a_1,b_1)
  recursive function SpAMM_tree_1d_dot_tree_1d_recur(a, b ) result(dot)

    type(SpAMM_tree_1d), pointer :: a,b

    real(SPAMM_KIND)             :: dot, dot0, dot1

    Dot=0

    if(.not.associated(a))return
    if(.not.associated(b))return

    if(a%frill%leaf)then

       dot = DOT_product( a%chunk(1:SBS), b%chunk(1:SBS) )

    else

       dot0=SpAMM_tree_1d_dot_tree_1d_recur( a%child_0, b%child_0 )
       dot1=SpAMM_tree_1d_dot_tree_1d_recur( a%child_1, b%child_1 )
       dot=dot0+dot1

    end if

  end function SpAMM_tree_1d_dot_tree_1d_recur

  !++NBODYTIMES:   SpAMM_init_random_tree_1d
  !++NBODYTIMES:     a => rand (wrapper)
  function SpAMM_random_tree_1d(M) result(randm)
    !
    integer,         intent(in)  :: M
    integer                      :: depth
    type(SpAMM_tree_1d), pointer :: randm
    real(SPAMM_KIND)             :: renorm

    randm => SpAMM_new_top_tree_1d(M)

    depth=0
    call init_random_seed()

    call SpAMM_random_unormalized_tree_1d_recur (randm, depth)

    ! normalize the vector ...

    renorm=SpAMM_one/sqrt(randm%frill%norm2)
    randm=>SpAMM_scalar_times_tree_1d(renorm, randm)

  end function SpAMM_random_tree_1d

  ! from the internet ...
  subroutine init_random_seed()

    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call RANDOM_seed(size = n)
    allocate(seed(n))

    call SYSTEM_clock(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call RANDOM_seed(PUT = seed)

    deallocate(seed)
  end subroutine init_random_seed

  !++NBODYTIMES:     SpAMM_random_unormalized_tree_1d_recur
  !++NBODYTIMES:       a_1 => rand (recursive)
  recursive subroutine SpAMM_random_unormalized_tree_1d_recur (randm, depth)

    type(SpAMM_tree_1d), pointer :: randm
    integer,          intent(in) :: depth
    integer                      :: hi,lo

    if(.not.associated(randm))return


    if(randm%frill%leaf)then

       randm%frill%needs_initialization = .false.

       lo=randm%frill%bndbx(0)
       hi=randm%frill%bndbx(1)

       !       randm%chunk(1:hi-lo+1)=SpAMM_one
       call RANDOM_number(randm%chunk(1:hi-lo+1))

    else

       ! child along [0]: [lo,mid] ...
       call SpAMM_random_unormalized_tree_1d_recur(SpAMM_construct_tree_1d_0(randm), depth+1 )
       ! child along [1]: [mid+1,hi] ...
       call SpAMM_random_unormalized_tree_1d_recur(SpAMM_construct_tree_1d_1(randm), depth+1 )

    end if

    ! merge & regarnish back up the tree ...
    call SpAMM_redecorate_tree_1d(randm)
    !
  end subroutine SpAMM_random_unormalized_tree_1d_recur

  !++NBODYTIMES:     SpAMM_scalar_times_tree_1d
  !++NBODYTIMES:       a_1 => alpha*a_1 wrapper)
  function SpAMM_scalar_times_tree_1d(alpha, a) result(d)

    type(SpAMM_tree_1d), pointer, intent(inout) :: a
    type(SpAMM_tree_1d), pointer                :: d
    real(SPAMM_KIND)                            :: alpha
    integer :: depth

    depth=0
    d=>a
    if(.not.associated(a))return

    call SpAMM_scalar_times_tree_1d_recur(alpha, d)

  end function SpAMM_scalar_times_tree_1d

  !++NBODYTIMES:     SpAMM_scalar_times_tree_1d_recur
  !++NBODYTIMES:       a_1 => alpha*a_1 (recursive)
  recursive subroutine SpAMM_scalar_times_tree_1d_recur(alpha, a)

    type(SpAMM_tree_1d), pointer :: a
    real(SPAMM_KIND)             :: alpha

    !    integer :: depth

    if(.not.associated(a))return

    if(a%frill%leaf)then

       a%frill%needs_initialization = .false.
       a%chunk=alpha*a%chunk
       a%frill%flops=a%frill%flops+SBS

    else
       ! child along [0]: [lo,mid] ...
       call SpAMM_scalar_times_tree_1d_recur(alpha, a%child_0 ) !,depth+1 )
       ! child along [1]: [mid+1,hi] ...
       call SpAMM_scalar_times_tree_1d_recur(alpha, a%child_1 ) !,depth+1 )
    end if

    ! merge & regarnish back up the tree ...
    call SpAMM_redecorate_tree_1d(a)
    !
  end subroutine SpAMM_scalar_times_tree_1d_recur

  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-ONE-D] ... [TREE-TWO-D X TREE-ONE-D] ...
  !++NBODYTIMES:     SpAMM_tree_2d_symm_times_tree_1d
  !++NBODYTIMES:     c_1 => alpha*c_1 + beta*(a_2.b_1) (wrapper)
  function SpAMM_tree_2d_symm_times_tree_1d(a, b, Tau, in_o) result(d)

    type(SpAMM_tree_2d_symm), pointer,           intent(IN)    :: A
    type(SpAMM_tree_1d),      pointer,           intent(IN)    :: B
    real(SPAMM_KIND),                            intent(IN)    :: Tau
    type(SpAMM_tree_1d),      pointer, optional                :: in_o
    type(SpAMM_tree_1d),      pointer                          :: D
    integer                                                    :: Depth
    real(SPAMM_KIND)                                           :: Tau2

    ! figure the starting conditions ...
    if(present(in_o))then
       d => in_o
    else
       d => null()
    end if
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared & modded SpAMM threshold, tau2 <- (tau*||A||*||B||)^2
    Tau2=Tau*Tau*a%frill%norm2*b%frill%norm2

    if(.not.associated(d))then
       ! instantiate a tree if no passed allocation
       d => SpAMM_new_top_tree_1d(a%frill%ndimn(1))
    end if

    ! set passed data for initialization
    call SpAMM_flip(d)

    Depth=0
    call SpAMM_tree_2d_symm_times_tree_1d_recur( d, A, B, Tau2, Depth )

    ! prune unused nodes ...
    call SpAMM_prune(d)

  end function SpAMM_tree_2d_symm_times_tree_1d

  recursive subroutine SpAMM_tree_2d_symm_times_tree_1d_recur( C, A, B, Tau2, Depth ) !<++NBODYTIMES|
    !                    c_1 => alpha*c_1 + beta*(aT_2.b_1) (recursive)                !<++NBODYTIMES|

    type(SpAMM_tree_2d_symm), pointer :: A !, INTENT(IN) :: A
    type(SpAMM_tree_1d),      pointer :: B !, INTENT(IN) :: B
    type(SpAMM_tree_1d),      pointer             :: C
    real(SPAMM_KIND),  intent(IN)                 :: Tau2
    integer                                       :: Depth
    type(SpAMM_tree_1d),      pointer             :: b0,b1
    type(SpAMM_tree_2d_symm), pointer             :: a00,a11,a01,a10

    logical :: tf

    if(c%frill%leaf )then ! Leaf condition ?
       if(c%frill%needs_initialization)then
          c%frill%needs_initialization = .false.
          c%chunk(1:SBS) = matmul( a%chunk(1:SBS,1:SBS), b%chunk(1:SBS) )
          c%frill%flops  = c%frill%flops + SBS2
       else
          c%chunk(1:SBS) = c%chunk(1:SBS) + matmul( a%chunk(1:SBS,1:SBS) , b%chunk(1:SBS) )
          c%frill%flops  = c%frill%flops + SBS2 + SBS
       end if
    else
       b0=>b%child_0;   b1=>b%child_1
       a00=>a%child_00; a11=>a%child_11
       a01=>a%child_01; a10=>a%child_10

       if(SpAMM_occlude( a00, b0, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_0(c), a00, b0, Tau2, Depth+1)
       if(SpAMM_occlude( a11, b1, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_1(c), a11, b1, Tau2, Depth+1)

       if(SpAMM_occlude( a01, b1, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_0(c), a01, b1, Tau2, Depth+1)
       if(SpAMM_occlude( a10, b0, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_1(c), a10, b0, Tau2, Depth+1)

    end if

    call SpAMM_redecorate_tree_1d(c)

  end subroutine SpAMM_tree_2d_symm_times_tree_1d_recur
  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-TWO-D] ... [TREE-TWO-D X TREE-TWO-D] ...


  !++NBODYTIMES:     SpAMM_scalar_times_tree_2d
  !++NBODYTIMES:       a_2 => alpha*a_2 wrapper)
  function SpAMM_scalar_times_tree_2d_symm(alpha, a) result(d)

    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: a
    type(SpAMM_tree_2d_symm), pointer                :: d
    real(SPAMM_KIND)                                 :: alpha
    integer :: depth

    d=>a
    if(.not.associated(a))return

    depth=0
    call SpAMM_scalar_times_tree_2d_symm_recur(alpha, d, depth)

  end function SpAMM_scalar_times_tree_2d_symm

  !++NBODYTIMES:     SpAMM_scalar_times_tree_2d_recur
  !++NBODYTIMES:       a_1 => alpha*a_1 (recursive)
  recursive subroutine SpAMM_scalar_times_tree_2d_symm_recur(alpha, a, depth)

    type(SpAMM_tree_2d_symm), pointer :: a
    real(SPAMM_KIND)                  :: alpha
    integer :: depth

    if(.not.associated(a))return

    if(a%frill%leaf)then

       a%frill%needs_initialization=.false.
       a%chunk(1:SBS,1:SBS)=alpha*a%chunk(1:SBS,1:SBS)
       a%frill%flops=a%frill%flops+SBS

    else

       ! child along [00]:
       call SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_00 ,depth+1 )
       ! child along [01]:
       call SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_01 ,depth+1 )
       ! child along [10]:
       call SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_10 ,depth+1 )
       ! child along [11]:
       call SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_11 ,depth+1 )

    end if

    ! merge & regarnish back up the tree ...
    call SpAMM_redecorate_tree_2d_symm(a)
    !
  end subroutine SpAMM_scalar_times_tree_2d_symm_recur


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (wrapper)
  function SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau, NT_O, In_O , stream_file_O) result(d)

    type(SpAMM_tree_2d_symm), pointer,           intent(IN)    :: A, B
    real(SPAMM_KIND),                            intent(IN)    :: Tau
    logical, optional,                           intent(IN)    :: NT_O
    type(SpAMM_tree_2d_symm), pointer, optional, intent(INOUT) :: In_O
    type(SpAMM_tree_2d_symm), pointer                          :: d
    integer                                                    :: Depth
    logical                                                    :: NT
    real(SPAMM_KIND)                                           :: Tau2
    character(LEN=*), optional     :: stream_file_O

#ifdef SpAMM_PRINT_STREAM

    integer :: maxi, maxj, maxk, i, j, k, Max_Depth

    type(SpAMM_cubes), pointer     :: SpAMM_Stream, Current
    real(SPAMM_KIND)               :: Opacity,  a_scale, b_scale, c_scale, abc_scale
    real(SPAMM_KIND)               :: MaxNorm,MinNorm,Emm,Bee


    real(SPAMM_KIND), dimension(:,:,:), allocatable :: Field

#endif

    ! figure the starting conditions ...
    if(present(in_O))then
       if(associated(in_o))then
          d => in_O
       else
          d => null()
       end if
    else
       d => null()
    end if

    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold
    Tau2=Tau*Tau

    if(present(NT_O))then
       NT=NT_O    ! If NT_O==FALSE, then A^t.B
    else
       NT=.true.  ! default is A.B
    end if

    if(.not.associated(d))then
       ! instantiate a tree if no passed allocation
       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn)
    end if

    ! set passed data for initialization
    call SpAMM_flip(d)

#ifdef SpAMM_PRINT_STREAM
    if(present(STREAM_FILE_O))then
       Build_Stream=.true.
       allocate(SpAMM_stream)
       Stream=>SpAMM_stream
       number_stream_elements = 0
       open(unit=45, file=trim(adjustl(stream_file_o))//'.norms', status='NEW')
    else
       Build_Stream=.false.
    end if
#endif

    Depth=0
    call SpAMM_tree_2d_symm_TIMES_tree_2d_symm_recur(d, A, B, Tau2, NT, Depth )

    ! prune unused nodes ...
    call SpAMM_prune(d)

#ifdef SpAMM_PRINT_STREAM
    if(.not.present(STREAM_FILE_O)) return

    ! Edit, this is mixed up in the plot file.  The easy fix is A->D, D->A.
    ! Then things make sense with visualization/plot-2.py

    write(45, "(A)") "Matrix A"
    call spamm_tree_print_leaves_2d_symm(D, file_unit=45)

    write(45, "(A)") "Matrix B"
    call spamm_tree_print_leaves_2d_symm(B, file_unit=45)

    write(45, "(A)") "Matrix C"
    call spamm_tree_print_leaves_2d_symm(A, file_unit=45)

    write(45, "(A)") "Product Space"
    write(45, "(I8)") SPAMM_CHUNK_SIZE
    do depth=0,64
       max_depth=depth
       if(SPAMM_CHUNK_SIZE*2**depth>=a%frill%NDimn(1))exit
    end do
    write(*,*)' max_depth = ', max_depth,SPAMM_CHUNK_SIZE* 2**max_depth

    !open(unit=44, file=trim(adjustl(stream_file_o))//'.vtk', status='NEW')

    !> write(44, "(A)") "# vtk DataFile Version 2.0"
    !> write(44, "(A)") "The product space"
    !> write(44, "(A)") "ASCII"
    !> write(44, "(A)") "DATASET POLYDATA"
    !> write(44, "(A,I30,A)") "POINTS ", number_stream_elements, " int"

    MaxNorm=-1D100
    MaxI=-100
    MaxJ=-100
    MaxK=-100
    MinNorm= 1D100
    Stream=>SpAMM_stream
    do while(associated(Stream%Next))
       i=Stream%Lw(1)+(Stream%Hi(1)-Stream%Lw(1)+1)/2
       j=Stream%Lw(2)+(Stream%Hi(2)-Stream%Lw(2)+1)/2
       k=Stream%Lw(3)+(Stream%Hi(3)-Stream%Lw(3)+1)/2
       !write(44, "(3I30)") i, j, k
       !> MaxI=MAX(MaxI,I)
       !> MaxJ=MAX(MaxJ,J)
       !> MaxK=MAX(MaxK,K)
       !> MaxNorm=MAX(MaxNorm,stream%size)
       !> MinNorm=MIN(MinNorm,stream%size)
       write(45, "(3ES20.10,3I8,ES20.10)") &
            stream%lw(1)+real(stream%hi(1)-stream%lw(1)+1)/2., &
            stream%lw(2)+real(stream%hi(2)-stream%lw(2)+1)/2., &
            stream%lw(3)+real(stream%hi(3)-stream%lw(3)+1)/2., &
            stream%hi(1)-stream%lw(1)+1, &
            stream%hi(2)-stream%lw(2)+1, &
            stream%hi(3)-stream%lw(3)+1, &
            log10(stream%size)
       Stream=>Stream%Next
    end do

    !close(44)
    close(45)

    !> WRITE(*,*)' MaxNorm = ',MinNorm, MaxNorm
    !> WRITE(*,*)' MaxIJK  = ',MaxI,MaxJ,MaxK
    !> WRITE(*,*)' MaxIJK  = ',MaxI*SPAMM_BLOCK_SIZE,MaxJ*SPAMM_BLOCK_SIZE,MaxK*SPAMM_BLOCK_SIZE

    !> ALLOCATE(FIELD(1:MaxI,1:MaxJ,1:MaxK))
    !> FIELD=0d0

    !> Stream=>SpAMM_stream
    !> DO WHILE(ASSOCIATED(Stream%Next))
    !>    i=ceiling(( Stream%Lw(1) + SpAMM_Half*( Stream%Hi(1)-Stream%Lw(1) ) )/SPAMM_BLOCK_SIZE )
    !>    j=ceiling(( Stream%Lw(2) + SpAMM_Half*( Stream%Hi(2)-Stream%Lw(2) ) )/SPAMM_BLOCK_SIZE )
    !>    k=ceiling(( Stream%Lw(3) + SpAMM_Half*( Stream%Hi(3)-Stream%Lw(3) ) )/SPAMM_BLOCK_SIZE )
    !>    Field(i,j,k)=stream%size
    !>    Stream=>Stream%Next
    !>END DO

    !> CALL VTK_write_scalar_3d(maxi,maxj,maxk,field,STREAM_FILE_O)
    !> DEALLOCATE(FIELD)

    Stream=>SpAMM_stream
    do while(associated(Stream%Next))
       Current=>Stream%Next
       deallocate(Stream)
       Stream=>Current
    end do
#endif

  end function SpAMM_tree_2d_symm_times_tree_2d_symm

  subroutine VTK_write_scalar_3d(ni,nj,nk, field, STREAM_FILE_O)

    integer, parameter          :: s=selected_real_kind(6)
    integer :: ni,nj,nk
    real(kind=SPAMM_KIND), intent(in), dimension(:,:,:) :: field
    character(len=*), optional                 :: STREAM_FILE_O
    character(len=1), parameter :: newline=achar(10)

    if(present(STREAM_FILE_O))then
       open(UNIT=44,FILE=trim(adjustl(STREAM_FILE_O))//'.vtk',STATUS='NEW')
    else
       stop ' Need to pass in file to open '
    end if

    write(44,'(A)')'# vtk DataFile Version 2.0'
    write(44,'(A)')'CT scan data of human heart, courtesy by Henk Mastenbroek RuG'
    write(44,'(A)')'ASCII'
    write(44,'(A)')' '
    write(44,'(A)')'DATASET STRUCTURED_POINTS'
    write(44,*)'DIMENSIONS',ni,nj,nk
    write(44,'(A)')'ORIGIN    1.000   1.000   1.000 '
    write(44,*)'SPACING', SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE, SPAMM_CHUNK_SIZE
    write(44,*)'POINT_DATA',ni*nj*nk
    write(44,'(A)')'SCALARS scalars float'
    write(44,'(A)')'LOOKUP_TABLE default'
    write(44,'(A)')' '
    write(44,*)real(field(1:ni,1:nj,1:nk),kind=s),newline
    close(unit=44)

  end subroutine VTK_write_scalar_3d


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
  !++NBODYTIMES:     c_2 => a_2 . b_2
  recursive subroutine SpAMM_tree_2d_symm_times_tree_2d_symm_recur( C, A, B, Tau2, NT, Depth )

    type(SpAMM_tree_2d_symm), pointer, intent(IN) :: A, B
    real(SPAMM_KIND),                  intent(IN) :: Tau2
    logical,                           intent(IN) :: NT
    integer,                           intent(IN) :: Depth
    type(SpAMM_tree_2d_symm), pointer             :: C
    type(SpAMM_tree_2d_symm), pointer             :: a00,a11,a01,a10
    type(SpAMM_tree_2d_symm), pointer             :: b00,b11,b01,b10
    type(SpAMM_tree_2d_symm), pointer             :: c00,c11,c01,c10

    if(c%frill%leaf) then ! Leaf condition ...

       if(c%frill%needs_initialization) then

          c%frill%needs_initialization = .false.

          if(NT)then
             c%chunk(1:SBS,1:SBS)=matmul(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))
          else
             c%chunk(1:SBS,1:SBS)=matmul(transpose(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
          end if

          c%frill%flops = SBS3

       else

          if(NT)then
             c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+matmul(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))
          else
             c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+matmul(transpose(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
          end if

          c%frill%flops = c%frill%flops + SBS2 + SBS3

       end if

#ifdef SpAMM_PRINT_STREAM
       if(Build_Stream)then
          if(NT)then
             Stream%Lw(1)=a%frill%bndbx(0,2)
             Stream%Lw(2)=a%frill%bndbx(0,1)
             Stream%Lw(3)=b%frill%bndbx(0,2)
             Stream%Hi(1)=a%frill%bndbx(1,2)
             Stream%Hi(2)=a%frill%bndbx(1,1)
             Stream%Hi(3)=b%frill%bndbx(1,2)
          else
             Stream%Lw(1)=a%frill%bndbx(0,1)
             Stream%Lw(2)=a%frill%bndbx(0,2)
             Stream%Lw(3)=b%frill%bndbx(0,2)
             Stream%Hi(1)=a%frill%bndbx(1,1)
             Stream%Hi(2)=a%frill%bndbx(1,2)
             Stream%Hi(3)=b%frill%bndbx(1,2)
          end if
          Stream%Levl=Depth
          Stream%Size=sqrt(a%frill%norm2*b%frill%norm2)
          number_stream_elements = number_stream_elements+1
          allocate(Stream%Next)
          Stream=>Stream%Next
          nullify(Stream%Next)
       end if
#endif

    else

       b00=>b%child_00; b11=>b%child_11; b01=>b%child_01; b10=>b%child_10
       a00=>a%child_00; a11=>a%child_11
       if(NT)then
          a01=>a%child_01; a10=>a%child_10
       else
          a01=>a%child_10; a10=>a%child_01
       end if

       ! first  pass, [m;0].[0;n]
       if(SpAMM_occlude( a00, b00, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_00(c),a00,b00,Tau2,NT,Depth+1)
       if(SpAMM_occlude( a10, b01, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_11(c),a10,b01,Tau2,NT,Depth+1)
       if(SpAMM_occlude( a00, b01, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_01(c),a00,b01,Tau2,NT,Depth+1)
       if(SpAMM_occlude( a10, b00, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_10(c),a10,b00,Tau2,NT,Depth+1)

       ! second  pass, [m;1].[1;n]
       if(SpAMM_occlude( a01, b10, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_00(c),a01,b10,Tau2,NT,Depth+1)
       if(SpAMM_occlude( a11, b11, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_11(c),a11,b11,Tau2,NT,Depth+1)
       if(SpAMM_occlude( a01, b11, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_01(c),a01,b11,Tau2,NT,Depth+1)
       if(SpAMM_occlude( a11, b10, Tau2 ) ) &
            call SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_10(c),a11,b10,Tau2,NT,Depth+1)
       !
    end if

    call SpAMM_redecorate_tree_2d_symm(c)

  end subroutine SpAMM_tree_2d_symm_times_tree_2d_symm_recur










!!$
!!$
!!$  FUNCTIOn SpAMM_tree_2d_symm_T_times_tree_2d_symm(a, b, Tau, alpha_O, beta_O, in_O) RESULT(d)
!!$
!!$    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
!!$    REAL(SPAMM_KIND),                            INTENT(IN)    :: Tau
!!$    REAL(SPAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
!!$    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: In_O
!!$    TYPE(SpAMM_tree_2d_symm), POINTER                          :: D
!!$    REAL(SPAMM_KIND)                                           :: alpha, beta
!!$    INTEGER                                                    :: Depth
!!$    REAL(SPAMM_KIND)                                           :: Tau2
!!$    ! figure the starting conditions ...
!!$    if(present(in_O))then
!!$       d => in_O
!!$    else
!!$       d => NULL()
!!$   end if
!!$    ! bail if we can ...
!!$    if(.not.associated(a))return
!!$    if(.not.associated(b))return
!!$
!!$    ! here is the squared threshold
!!$    Tau2=Tau*Tau
!!$
!!$    ! need a new tree? then instantiate one ...
!!$    if(.not.associated(d))&
!!$       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn)
!!$
!!$    if(present(alpha_O))then
!!$       d => SpAMM_scalar_times_tree_2d_symm(alpha_O, d)
!!$   end if
!!$
!!$    beta =SpAMM_one
!!$    if(present( beta_O))beta = beta_O
!!$
!!$    CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(d, A, B, Tau2, Depth, beta )
!!$
!!$  END FUNCTION SpAMM_tree_2d_symm_T_times_tree_2d_symm
!!$
!!$
!!$  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
!!$  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (recursive)
!!$  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(C, A, B, Tau2, Depth,  beta )
!!$
!!$    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
!!$    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
!!$    REAL(SPAMM_KIND),                  INTENT(IN)    :: beta
!!$    REAL(SPAMM_KIND)                                 :: Tau2
!!$    INTEGER                                          :: Depth
!!$    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c10,c11
!!$
!!$    integer, dimension(1:2) :: ahi,alo , bhi,blo , chi , clo, athi , atlo
!!$
!!$    if(.not.associated(a))return
!!$    if(.not.associated(b))return
!!$
!!$    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
!!$    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return
!!$
!!$    IF(c%frill%leaf )THEN ! Leaf condition ?
!!$
!!$!       WRITE(*,*)'before ',SQRT(SUM(c%chunk**2)),' beta = ',beta,SQRT(SUM(a%chunk**2)),SQRT(SUM(b%chunk**2))
!!$
!!$ !      blo=b%frill%bndbx(0,:)
!!$ !      bhi=b%frill%bndbx(1,:)
!!$ !      clo=c%frill%bndbx(0,:)
!!$ !      chi=c%frill%bndbx(1,:)
!!$ !      alo=a%frill%bndbx(0,:)
!!$ !      ahi=a%frill%bndbx(1,:)
!!$
!!$ !       aTlo(1)=a%frill%bndbx(0,2)
!!$ !      aTlo(2)=a%frill%bndbx(0,1)
!!$ !      aThi(1)=a%frill%bndbx(1,2)
!!$ !      aThi(2)=a%frill%bndbx(1,1)
!!$
!!$
!!$!       WRITE(*,67)clo(1),chi(1), clo(2),chi(2), &
!!$!                  aTlo(1),aThi(1), aTlo(2),aThi(2), &
!!$!                  blo(1),bhi(1), blo(2),bhi(2)
!!$
!!$!       WRITE(*,77)chi(1)-clo(1)+1, chi(2)-clo(2)+1, &
!!$!                  ahi(1)-alo(1)+1, ahi(2)-alo(2)+1, &
!!$!                  bhi(1)-blo(1)+1, bhi(2)-blo(2)+1
!!$
!!$67     format('[',I3,'-',I3,', ',I3,'-',I3,'] = [',I3,'-',I3,', ',I3,'-',I3,']^t x [',I3,'-',I3,', ',I3,'-',I3,']')
!!$77     format('[',I3,', ',I3,'] = [',I3,', ',I3,'] x [',I3,', ',I3,']')
!!$
!!$!       c%chunk(1:(chi(1)-clo(1)+1), 1:(chi(2)-clo(2)+1)) = &
!!$!       c%chunk(1:(chi(1)-clo(1)+1), 1:(chi(2)-clo(2)+1)) + &
!!$!       beta*MATMUL( TRANSPOSE( a%chunk( 1:ahi(1)-alo(1)+1 , 1:ahi(2)-alo(2)+1 ) ), &
!!$!                   b%chunk(1:(bhi(1)-blo(1)+1), 1:(bhi(2)-blo(2)+1)))
!!$
!!$       c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+ &
!!$               beta*MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
!!$
!!$!       WRITE(*,*)'after ',SQRT(SUM(c%chunk**2))
!!$
!!$       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2
!!$
!!$   ELSE
!!$
!!$       ! find some memory ...
!!$       c00=>SpAMM_construct_tree_2d_symm_00(c)
!!$       c01=>SpAMM_construct_tree_2d_symm_01(c)
!!$       c10=>SpAMM_construct_tree_2d_symm_10(c)
!!$       c11=>SpAMM_construct_tree_2d_symm_11(c)
!!$
!!$       ! a first pass ...
!!$!       WRITE(*,*)' 00 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c00, a%child_00, b%child_00, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$!       WRITE(*,*)' 01=00:01 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c01, a%child_00, b%child_01, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$!       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c10, a%child_10, b%child_00, &
!!$!                                             Tau2, Depth+1, beta )
!!$
!!$ !      WRITE(*,*)' 10=01^t:00 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c10, a%child_01, b%child_00, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$!       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c11, a%child_10, b%child_01, &
!!$!                                             Tau2, Depth+1, beta )
!!$
!!$!       WRITE(*,*)' 11=01^t:01 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c11, a%child_01, b%child_01, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$!       WRITE(*,*)' ------------------'
!!$
!!$       ! ... & another pass
!!$!       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c00, a%child_01, b%child_10, &
!!$!                                             Tau2, Depth+1, beta )
!!$!       WRITE(*,*)' 00=10^t:10 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c00, a%child_10, b%child_10, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$
!!$!       WRITE(*,*)' 01=10^t:11 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c01, a%child_10, b%child_11, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$!       WRITE(*,*)' 10=11:10 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c10, a%child_11, b%child_10, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$!       WRITE(*,*)' 11=11:11 '
!!$       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c11, a%child_11, b%child_11, &
!!$                                             Tau2, Depth+1, beta )
!!$
!!$       !
!!$   END IF
!!$
!!$    CALL SpAMM_redecorate_tree_2d_symm(c)
!!$
!!$  END SUBROUTINE SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur



end module spamm_nbdyalgbra_times
