
#define SpAMM_PRINT_STREAM
module spamm_nbdyalgbra_times 

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals

  implicit none

  TYPE SpAMM_cubes
     TYPE(SpAMM_cubes), POINTER  :: Next
     INTEGER                     :: Levl
     REAL(SpAMM_Kind)            :: Size
     INTEGER,     DIMENSION(3)   :: Lw, Hi
     TYPE(SpAMM_tree_2d_symm),POINTER    :: A,B,C
  END TYPE SpAMM_cubes

  TYPE(SpAMM_cubes), POINTER     :: Stream


CONTAINS



  !++NBODYTIMES: SpAMM generalized n-body algebras for times ____ NBODYTIMES _________________
  !++NBODYTIMES: generalized products, dots, contractions & convolutions (X)
  !++NBODYTIMES:   ... [TREE-ONE-D X TREE-ONE-D] ... [TREE-ONE-D X TREE-ONE-D] ...   
  !++NBODYTIMES:   SpAMM_tree_1d_dot_tree_1d_recur
  !++NBODYTIMES:     dot = (a_1,b_1) 
  RECURSIVE FUNCTION SpAMM_tree_1d_dot_tree_1d_recur(a, b ) result(dot)

    TYPE(SpAMM_tree_1d), POINTER :: a,b

    REAL(SpAMM_KIND)             :: dot, dot0, dot1

    Dot=SpAMM_Zero

    if(.not.associated(a))return
    if(.not.associated(b))return

    if(a%frill%leaf)then


       dot = DOT_PRODUCT( a%chunk(1:SBS), b%chunk(1:SBS) )

    else

       dot0=SpAMM_tree_1d_dot_tree_1d_recur( a%child_0, b%child_0 )
       dot1=SpAMM_tree_1d_dot_tree_1d_recur( a%child_1, b%child_1 )
       dot=dot0+dot1

    endif

  END FUNCTION SpAMM_tree_1d_dot_tree_1d_recur

  !++NBODYTIMES:   SpAMM_init_random_tree_1d
  !++NBODYTIMES:     a => rand (wrapper)
  function SpAMM_random_tree_1d(M) result (randm)
    !
    integer,         intent(in)  :: M
    integer                      :: depth
    type(SpAMM_tree_1d), pointer :: randm
    real(SpAMM_KIND)             :: renorm

    randm => SpAMM_new_top_tree_1d(M)

    depth=0
    CALL init_random_seed()

    CALL SpAMM_random_unormalized_tree_1d_recur (randm, depth)

    ! normalize the vector ...

    renorm=SpAMM_one/sqrt(randm%frill%norm2)
    randm=>SpAMM_scalar_times_tree_1d(renorm, randm)

  end function SpAMM_random_tree_1d

  ! from the internet ...
  subroutine init_random_seed()
    
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    CALL SYSTEM_CLOCK(COUNT=clock)
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    
    DEALLOCATE(seed)
  end subroutine init_random_seed

  !++NBODYTIMES:     SpAMM_random_unormalized_tree_1d_recur
  !++NBODYTIMES:       a_1 => rand (recursive)
  recursive subroutine SpAMM_random_unormalized_tree_1d_recur (randm, depth)

    type(SpAMM_tree_1d), pointer :: randm
    integer,          intent(in) :: depth 
    INTEGER                      :: hi,lo

    if(.not.associated(randm))return

    
    IF(randm%frill%leaf)THEN

       randm%frill%init = .FALSE.
       
       lo=randm%frill%bndbx(0)
       hi=randm%frill%bndbx(1)

!       randm%chunk(1:hi-lo+1)=SpAMM_one
       CALL RANDOM_NUMBER(randm%chunk(1:hi-lo+1))

    ELSE

       ! child along [0]: [lo,mid] ... 
       CALL SpAMM_random_unormalized_tree_1d_recur(SpAMM_construct_tree_1d_0(randm), depth+1 )
       ! child along [1]: [mid+1,hi] ... 
       CALL SpAMM_random_unormalized_tree_1d_recur(SpAMM_construct_tree_1d_1(randm), depth+1 )

    ENDIF
    
    ! merge & regarnish back up the tree ...
    CALL SpAMM_redecorate_tree_1d(randm)
    !
  end subroutine SpAMM_random_unormalized_tree_1d_recur

  !++NBODYTIMES:     SpAMM_scalar_times_tree_1d
  !++NBODYTIMES:       a_1 => alpha*a_1 wrapper)
  FUNCTION SpAMM_scalar_times_tree_1d(alpha, a) RESULT(d)
 
    type(SpAMM_tree_1d), pointer, intent(inout) :: a
    type(SpAMM_tree_1d), pointer                :: d
    real(SpAMM_KIND)                            :: alpha
    integer :: depth

    depth=0
    d=>a
    if(.not.associated(a))return

    CALL SpAMM_scalar_times_tree_1d_recur(alpha, d)

  END FUNCTION SpAMM_scalar_times_tree_1d

  !++NBODYTIMES:     SpAMM_scalar_times_tree_1d_recur
  !++NBODYTIMES:       a_1 => alpha*a_1 (recursive)
  recursive subroutine SpAMM_scalar_times_tree_1d_recur(alpha, a)

    type(SpAMM_tree_1d), pointer :: a
    real(SpAMM_KIND)             :: alpha

!    integer :: depth

    if(.not.associated(a))return

    IF(a%frill%leaf)THEN

       a%frill%init = .FALSE.
       a%chunk=alpha*a%chunk
       a%frill%flops=a%frill%flops+SBS

    ELSE
       ! child along [0]: [lo,mid] ... 
       CALL SpAMM_scalar_times_tree_1d_recur(alpha, a%child_0 ) !,depth+1 )
       ! child along [1]: [mid+1,hi] ... 
       CALL SpAMM_scalar_times_tree_1d_recur(alpha, a%child_1 ) !,depth+1 )
    ENDIF  

    ! merge & regarnish back up the tree ...
    CALL SpAMM_redecorate_tree_1d(a)
    !
  end subroutine SpAMM_scalar_times_tree_1d_recur

  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-ONE-D] ... [TREE-TWO-D X TREE-ONE-D] ...   
  !++NBODYTIMES:     SpAMM_tree_2d_symm_times_tree_1d
  !++NBODYTIMES:     c_1 => alpha*c_1 + beta*(a_2.b_1) (wrapper)
  FUNCTION SpAMM_tree_2d_symm_times_tree_1d(a, b, Tau, in_o) RESULT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A
    TYPE(SpAMM_tree_1d),      POINTER,           INTENT(IN)    :: B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    TYPE(SpAMM_tree_1d),      POINTER, OPTIONAL  :: in_o
    TYPE(SpAMM_tree_1d),      POINTER                          :: D
    INTEGER                                                    :: Depth
    REAL(SpAMM_KIND)                                           :: Tau2

    ! figure the starting conditions ...
    if(present(in_o))then
       d => in_o
    else
       d => NULL()
    endif
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    if(.not.associated(d))then
       ! instantiate a tree if no passed allocation
       d => SpAMM_new_top_tree_1d(a%frill%ndimn(1))
    endif

    ! set passed data for initialization
    CALL SpAMM_flip(d)

    Depth=0
    CALL SpAMM_tree_2d_symm_times_tree_1d_recur( d, A, B, Tau2, Depth )

    ! prune unused nodes ... 
    CALL SpAMM_prune(d)

  END FUNCTION SpAMM_tree_2d_symm_times_tree_1d

  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_1d_recur( C, A, B, Tau2, Depth ) !<++NBODYTIMES|
   !                    c_1 => alpha*c_1 + beta*(aT_2.b_1) (recursive)                !<++NBODYTIMES|    

    TYPE(SpAMM_tree_2d_symm), POINTER :: A !, INTENT(IN) :: A
    TYPE(SpAMM_tree_1d),      POINTER :: B !, INTENT(IN) :: B
    TYPE(SpAMM_tree_1d),      POINTER             :: C
    REAL(SpAMM_KIND),  INTENT(IN)                 :: Tau2
    INTEGER                                       :: Depth
    TYPE(SpAMM_tree_1d),      POINTER             :: b0,b1
    TYPE(SpAMM_tree_2d_symm), POINTER             :: a00,a11,a01,a10

    logical :: tf

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       IF( c%frill%init )THEN

          c%frill%init   = .FALSE.
          c%chunk(1:SBS) = MATMUL( a%chunk(1:SBS,1:SBS), b%chunk(1:SBS) )
          c%frill%flops  = c%frill%flops + SBS2
          
       ELSE

          c%chunk(1:SBS) = c%chunk(1:SBS) + MATMUL( a%chunk(1:SBS,1:SBS) , b%chunk(1:SBS) )
          c%frill%flops  = c%frill%flops + SBS2 + SBS
          
       ENDIF

    ELSE

        b0=>b%child_0;   b1=>b%child_1
       a00=>a%child_00; a11=>a%child_11 
       a01=>a%child_01; a10=>a%child_10

       IF( SpAMM_occlude( a00, b0, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_0(c), a00, b0, Tau2, Depth+1)
       IF( SpAMM_occlude( a11, b1, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_1(c), a11, b1, Tau2, Depth+1)

       IF( SpAMM_occlude( a01, b1, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_0(c), a01, b1, Tau2, Depth+1)
       IF( SpAMM_occlude( a10, b0, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_1d_recur(SpAMM_construct_tree_1d_1(c), a10, b0, Tau2, Depth+1)

    ENDIF

    CALL SpAMM_redecorate_tree_1d(c)

  END SUBROUTINE SpAMM_tree_2d_symm_times_tree_1d_recur
  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-TWO-D] ... [TREE-TWO-D X TREE-TWO-D] ...   


  !++NBODYTIMES:     SpAMM_scalar_times_tree_2d
  !++NBODYTIMES:       a_2 => alpha*a_2 wrapper)
  FUNCTION SpAMM_scalar_times_tree_2d_symm(alpha, a) RESULT(d)
 
    type(SpAMM_tree_2d_symm), pointer, intent(inout) :: a
    type(SpAMM_tree_2d_symm), pointer                :: d
    real(SpAMM_KIND)                                 :: alpha
    integer :: depth

    d=>a
    if(.not.associated(a))return

    depth=0
    CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, d, depth)

  END FUNCTION SpAMM_scalar_times_tree_2d_symm

  !++NBODYTIMES:     SpAMM_scalar_times_tree_2d_recur
  !++NBODYTIMES:       a_1 => alpha*a_1 (recursive)
  recursive subroutine SpAMM_scalar_times_tree_2d_symm_recur(alpha, a, depth)

    type(SpAMM_tree_2d_symm), pointer :: a
    real(SpAMM_KIND)                  :: alpha
    integer :: depth

    if(.not.associated(a))return

    IF(a%frill%leaf)THEN

       a%frill%init=.FALSE.
       a%chunk(1:SBS,1:SBS)=alpha*a%chunk(1:SBS,1:SBS)
       a%frill%flops=a%frill%flops+SBS

    ELSE

       ! child along [00]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_00 ,depth+1 )
       ! child along [01]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_01 ,depth+1 )
       ! child along [10]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_10 ,depth+1 )
       ! child along [11]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_11 ,depth+1 )

    ENDIF  

    ! merge & regarnish back up the tree ...
    CALL SpAMM_redecorate_tree_2d_symm(a)
    !
  end subroutine SpAMM_scalar_times_tree_2d_symm_recur


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (wrapper)
#ifdef SpAMM_PRINT_STREAM
  FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau, NT_O, In_O , stream_file_O) RESULT(d)
    !
#else
  FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau, NT_O, In_O ) RESULT(d)
#endif
    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    LOGICAL, OPTIONAL,                           INTENT(IN)    :: NT_O
    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: In_O
    TYPE(SpAMM_tree_2d_symm), POINTER                          :: d
    INTEGER                                                    :: Depth
    LOGICAL                                                    :: NT
    REAL(SpAMM_KIND)                                           :: Tau2
    
#ifdef SpAMM_PRINT_STREAM

    integer :: maxi, maxj, maxk, i, j, k, Max_Depth

    TYPE(SpAMM_cubes), POINTER     :: SpAMM_Stream
    REAL(SpAMM_kind)               :: Opacity,  a_scale, b_scale, c_scale, abc_scale
    REAL(SpAMM_kind)               :: MaxNorm,MinNorm,Emm,Bee

    CHARACTER(LEN=*), OPTIONAL     :: stream_file_O
    
    REAL(SpAMM_Kind), DIMENSION(:,:,:), ALLOCATABLE :: Field

#endif

    ! figure the starting conditions ...
    if(present(in_O))then       
       if(associated(in_o))then
          d => in_O
       else
          d => NULL()
       endif
    else
       d => NULL()
    endif

    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    if(present(NT_O))then
       NT=NT_O    ! If NT_O==FALSE, then A^t.B
    else
       NT=.TRUE.  ! default is A.B
    endif

    if(.not.associated(d))then
       ! instantiate a tree if no passed allocation
       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn)
    endif
 
    ! set passed data for initialization   
    CALL SpAMM_flip(d)
    
#ifdef SpAMM_PRINT_STREAM
    ALLOCATE(SpAMM_stream)
    Stream=>SpAMM_stream
#endif

    Depth=0
    CALL SpAMM_tree_2d_symm_TIMES_tree_2d_symm_recur(d, A, B, Tau2, NT, Depth )

    ! prune unused nodes ... 
    CALL SpAMM_prune(d)

#ifdef SpAMM_PRINT_STREAM
    IF(.NOT.PRESENT(STREAM_FILE_O))RETURN

    do depth=0,64       
       max_depth=depth
       if(SPAMM_BLOCK_SIZE*2**depth>=a%frill%NDimn(1))exit
    enddo
    WRITE(*,*)' max_depth = ', max_depth,SPAMM_BLOCK_SIZE* 2**max_depth

    MaxNorm=-1D100
    MaxI=-100
    MaxJ=-100
    MaxK=-100
    MinNorm= 1D100
    Stream=>SpAMM_stream
    DO WHILE(ASSOCIATED(Stream%Next))
       IF(Stream%Levl==max_depth)THEN
          MaxNorm=MAX(MaxNorm,SQRT(stream%a%frill%norm2))
          MaxNorm=MAX(MaxNorm,SQRT(stream%b%frill%norm2))
          MaxNorm=MAX(MaxNorm,SQRT(stream%c%frill%norm2))
          MinNorm=MIN(MinNorm,SQRT(stream%a%frill%norm2))
          MinNorm=MIN(MinNorm,SQRT(stream%b%frill%norm2))
          MinNorm=MIN(MinNorm,SQRT(stream%c%frill%norm2))

          i=ceiling( ( Stream%Lw(1) + SpAMM_Half*( Stream%Hi(1)-Stream%Lw(1) ) )/SPAMM_BLOCK_SIZE )
          j=ceiling( ( Stream%Lw(2) + SpAMM_Half*( Stream%Hi(2)-Stream%Lw(2) ) )/SPAMM_BLOCK_SIZE )
          k=ceiling( ( Stream%Lw(3) + SpAMM_Half*( Stream%Hi(3)-Stream%Lw(3) ) )/SPAMM_BLOCK_SIZE )
          MaxI=MAX(MaxI,I)
          MaxJ=MAX(MaxJ,J)
          MaxK=MAX(MaxK,K)
       ENDIF
       Stream=>Stream%Next 
    ENDDO

    WRITE(*,*)' MaxNorm = ',MinNorm, MaxNorm
    WRITE(*,*)' MaxIJK  = ',MaxI,MaxJ,MaxK
    WRITE(*,*)' MaxIJK  = ',MaxI*SPAMM_BLOCK_SIZE,MaxJ*SPAMM_BLOCK_SIZE,MaxK*SPAMM_BLOCK_SIZE

    ALLOCATE(FIELD(1:MaxI,1:MaxJ,1:MaxK))
    FIELD=0d0


!!$
!!$
!!$    MaxNorm=LOG10(MaxNorm)
!!$    MinNorm=LOG10(MinNorm)
!!$    WRITE(*,*)' MaxNorm = ',MinNorm, MaxNorm
!!$
!!$    Emm=1D0/(MaxNorm-MinNorm)
!!$    Bee=Emm*MinNorm
!!$    WRITE(*,*)' Bee, Emm = ',Bee,Emm
!!$
!!$
!!$    WRITE(44,*)'Graphics3D[{EdgeForm[],'

    Stream=>SpAMM_stream
    DO WHILE(ASSOCIATED(Stream%Next))
       IF(Stream%Levl==max_depth)THEN
          Opacity=1D0
!       ELSE
!          Opacity=0.4D0*DBLE(Stream%Levl+0.1d0)/DBLE(max_depth)
!       ENDIF

       IF(ASSOCIATED(Stream%Next%next))THEN

!!$          a_scale  =-Bee+Emm*LOG10(SQRT(stream%a%frill%norm2))
!!$          b_scale  =-Bee+Emm*LOG10(SQRT(stream%b%frill%norm2))
!!$          c_scale  =-Bee+Emm*LOG10(SQRT(stream%c%frill%norm2))
!!$          abc_scale=-Bee+Emm*LOG10(stream%size)

!!$          WRITE(44,111)a_scale,transpose(stream%a%frill%bndbx),  &
!!$                       b_scale,transpose(stream%b%frill%bndbx),  &
!!$                       c_scale,transpose(stream%c%frill%bndbx),  &
!!$                       abc_scale,Stream%Lw(1),Stream%Lw(2),Stream%Lw(3),Stream%Hi(1),Stream%Hi(2),Stream%Hi(3)
!!$111 FORMAT("Opacity[",F10.5,"], Cuboid[{",I6,",",I6,", -5.0 },{ ",I6,",",I6,", -5.01 }] , ",  &
!!$           "Opacity[",F10.5,"], Cuboid[{ -5.0 ,",I6,",",I6,"},{ -5.01 ,",I6,",",I6, "}] , ",  &
!!$           "Opacity[",F10.5,"], Cuboid[{",I6,", -5.0 ,",I6,"},{",I6,", -5.01 ,",I6, "}] , ",  &
!!$           "Opacity[",F10.5,"], Cuboid[{",I6,",",I6,",",I6,"},{",I6,",",I6,",",I6,"}] , ")


!!$          WRITE(44,111)a_scale,dble(transpose(stream%a%frill%bndbx)),  &
!!$                       b_scale,dble(transpose(stream%b%frill%bndbx)),  &
!!$                       c_scale,dble(transpose(stream%c%frill%bndbx)),  &
!!$                       abc_scale,Opacity,dble(Stream%Lw(1:3)),dble(Stream%Hi(1:3))
!!$
!!$111 FORMAT('ColorData["Pastel"][',F10.5,"], Cuboid[{",F10.5,",",F10.5,", -5.0 },{ ",F10.5,",",F10.5,", -5.01 }] , ",  &
!!$           'ColorData["Pastel"][',F10.5,"], Cuboid[{ -5.0 ,",F10.5,",",F10.5,"},{ -5.01 ,",F10.5,",",F10.5, "}] , ",  &
!!$           'ColorData["Pastel"][',F10.5,"], Cuboid[{",F10.5,", -5.0 ,",F10.5,"},{",F10.5,", -5.01 ,",F10.5, "}] , ",  &
!!$           'ColorData["Pastel"][',F10.5,"], Opacity[",F10.5,"], Cuboid[{",F10.5,",",F10.5,",",F10.5,"},{",F10.5,",",F10.5,",",F10.5,"}] , ")

          i=ceiling(( Stream%Lw(1) + SpAMM_Half*( Stream%Hi(1)-Stream%Lw(1) ) )/SPAMM_BLOCK_SIZE )
          j=ceiling(( Stream%Lw(2) + SpAMM_Half*( Stream%Hi(2)-Stream%Lw(2) ) )/SPAMM_BLOCK_SIZE )
          k=ceiling(( Stream%Lw(3) + SpAMM_Half*( Stream%Hi(3)-Stream%Lw(3) ) )/SPAMM_BLOCK_SIZE )

          Field(i,j,k)=stream%size

!          WRITE(44,111)abc_scale,dble(Stream%Lw(1:3)),dble(Stream%Hi(1:3))

111 FORMAT('ColorData["Pastel"][',F10.5,"], Cuboid[{",F10.5,",",F10.5,",",F10.5,"},{",F10.5,",",F10.5,",",F10.5,"}] , ")


!!$          WRITE(44,111)a_scale,transpose(stream%a%frill%bndbx),  &
!!$                       b_scale,transpose(stream%b%frill%bndbx),  &
!!$                       c_scale,transpose(stream%c%frill%bndbx),  &
!!$                       abc_scale,Opacity,Stream%Lw(1),Stream%Lw(2),Stream%Lw(3),Stream%Hi(1),Stream%Hi(2),Stream%Hi(3)
!!$111 FORMAT('ColorData["Pastel"][',F10.5,"], Cuboid[{",I6,",",I6,", -5.0 },{ ",I6,",",I6,", -5.01 }] , ",  &
!!$           'ColorData["Pastel"][',F10.5,"], Cuboid[{ -5.0 ,",I6,",",I6,"},{ -5.01 ,",I6,",",I6, "}] , ",  &
!!$           'ColorData["Pastel"][',F10.5,"], Cuboid[{",I6,", -5.0 ,",I6,"},{",I6,", -5.01 ,",I6, "}] , ",  &
!!$           'ColorData["Pastel"][',F10.5,"], Opacity[",F10.5,"], Cuboid[{",I6,",",I6,",",I6,"},{",I6,",",I6,",",I6,"}] , ")
!!$



       ENDIF


    ENDIF
 
       Stream=>Stream%Next 
    ENDDO


    CALL VTK_write_scalar_3d(maxi,maxj,maxk,field,STREAM_FILE_O)
    DEALLOCATE(FIELD)


!    Stream=>SpAMM_stream
!    DO WHILE(ASSOCIATED(Stream%Next))       
!    ENDDO
#endif 
 
  END FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm

  subroutine VTK_write_scalar_3d(ni,nj,nk, field, STREAM_FILE_O)

    integer, parameter          :: s=selected_real_kind(6)
    integer :: ni,nj,nk
    real(kind=SpAMM_Kind), intent(in), dimension(:,:,:) :: field
    character(len=*), optional                 :: STREAM_FILE_O
    character(len=1), parameter :: newline=achar(10)

    IF(PRESENT(STREAM_FILE_O))THEN
       OPEN(UNIT=44,FILE=TRIM(ADJUSTL(STREAM_FILE_O))//'.vtk',STATUS='NEW')
    ELSE
       STOP ' Need to pass in file to open '
    ENDIF

    WRITE(44,'(A)')'# vtk DataFile Version 2.0'
    WRITE(44,'(A)')'CT scan data of human heart, courtesy by Henk Mastenbroek RuG'
    WRITE(44,'(A)')'ASCII'
    WRITE(44,'(A)')' '
    WRITE(44,'(A)')'DATASET STRUCTURED_POINTS'
    WRITE(44,*)'DIMENSIONS',ni,nj,nk
    WRITE(44,'(A)')'ORIGIN    1.000   1.000   1.000 '
    WRITE(44,*)'SPACING', SPAMM_BLOCK_SIZE,SPAMM_BLOCK_SIZE,SPAMM_BLOCK_SIZE
    WRITE(44,*)'POINT_DATA',ni*nj*nk
    WRITE(44,'(A)')'SCALARS scalars float'
    WRITE(44,'(A)')'LOOKUP_TABLE default'
    WRITE(44,'(A)')' '
    write(44,*)real(field(1:ni,1:nj,1:nk),kind=s),newline
    CLOSE(unit=44)

  end subroutine VTK_write_scalar_3d



  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
  !++NBODYTIMES:     c_2 => a_2 . b_2 
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur( C, A, B, Tau2, NT, Depth )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN) :: A, B
    REAL(SpAMM_KIND),                  INTENT(IN) :: Tau2
    LOGICAL,                           INTENT(IN) :: NT
    INTEGER,                           INTENT(IN) :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER             :: C
    TYPE(SpAMM_tree_2d_symm), POINTER             :: a00,a11,a01,a10
    TYPE(SpAMM_tree_2d_symm), POINTER             :: b00,b11,b01,b10
    TYPE(SpAMM_tree_2d_symm), POINTER             :: c00,c11,c01,c10

    IF( c%frill%leaf )THEN ! Leaf condition ... 

#ifdef SpAMM_PRINT_STREAM
       Stream%a=>a
       Stream%b=>b
       Stream%c=>c
#endif

       IF( c%frill%init )THEN

          c%frill%init = .FALSE.

          IF(NT)THEN
             c%chunk(1:SBS,1:SBS)=MATMUL(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))
          ELSE
             c%chunk(1:SBS,1:SBS)=MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
          ENDIF

          c%frill%flops = SBS3

       ELSE

          IF(NT)THEN
             c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+MATMUL(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))
          ELSE
             c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
          ENDIF

          c%frill%flops = c%frill%flops + SBS2 + SBS3

       ENDIF

#ifdef SpAMM_PRINT_STREAM


       IF(NT)THEN
          Stream%Lw(1)=a%frill%bndbx(0,2)
          Stream%Lw(2)=a%frill%bndbx(0,1)
          Stream%Lw(3)=b%frill%bndbx(0,2)
          Stream%Hi(1)=a%frill%bndbx(1,2)
          Stream%Hi(2)=a%frill%bndbx(1,1)
          Stream%Hi(3)=b%frill%bndbx(1,2)
       ELSE
          Stream%Lw(1)=a%frill%bndbx(0,1)
          Stream%Lw(2)=a%frill%bndbx(0,2)
          Stream%Lw(3)=b%frill%bndbx(0,2)
          Stream%Hi(1)=a%frill%bndbx(1,1)
          Stream%Hi(2)=a%frill%bndbx(1,2)
          Stream%Hi(3)=b%frill%bndbx(1,2)
       ENDIF

       Stream%Levl=Depth
       Stream%Size=SQRT(a%frill%norm2*b%frill%norm2)

       ALLOCATE(Stream%Next)
       Stream=>Stream%Next
       NULLIFY(Stream%Next)
#endif


   ELSE


       b00=>b%child_00; b11=>b%child_11; b01=>b%child_01; b10=>b%child_10
       a00=>a%child_00; a11=>a%child_11
       IF(NT)THEN  
          a01=>a%child_01; a10=>a%child_10
       ELSE
          a01=>a%child_10; a10=>a%child_01
       ENDIF




       ! first  pass, [m;0].[0;n]
       IF( SpAMM_occlude( a00, b00, Tau2 ) ) & 
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_00(c),a00,b00,Tau2,NT,Depth+1)
       IF( SpAMM_occlude( a10, b01, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_11(c),a10,b01,Tau2,NT,Depth+1)
       IF( SpAMM_occlude( a00, b01, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_01(c),a00,b01,Tau2,NT,Depth+1)
       IF( SpAMM_occlude( a10, b00, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_10(c),a10,b00,Tau2,NT,Depth+1)

       ! second  pass, [m;1].[1;n]
       IF( SpAMM_occlude( a01, b10, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_00(c),a01,b10,Tau2,NT,Depth+1)
       IF( SpAMM_occlude( a11, b11, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_11(c),a11,b11,Tau2,NT,Depth+1)
       IF( SpAMM_occlude( a01, b11, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_01(c),a01,b11,Tau2,NT,Depth+1)
       IF( SpAMM_occlude( a11, b10, Tau2 ) ) &
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(SpAMM_construct_tree_2d_symm_10(c),a11,b10,Tau2,NT,Depth+1)
       !
    ENDIF

    CALL SpAMM_redecorate_tree_2d_symm(c)

  END SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur










!!$
!!$
!!$  FUNCTIOn SpAMM_tree_2d_symm_T_times_tree_2d_symm(a, b, Tau, alpha_O, beta_O, in_O) RESULT(d)
!!$
!!$    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
!!$    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
!!$    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
!!$    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: In_O
!!$    TYPE(SpAMM_tree_2d_symm), POINTER                          :: D
!!$    REAL(SpAMM_KIND)                                           :: alpha, beta
!!$    INTEGER                                                    :: Depth
!!$    REAL(SpAMM_KIND)                                           :: Tau2
!!$    ! figure the starting conditions ...
!!$    if(present(in_O))then
!!$       d => in_O
!!$    else
!!$       d => NULL()
!!$    endif
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
!!$    endif
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
!!$    REAL(SpAMM_KIND),                  INTENT(IN)    :: beta
!!$    REAL(SpAMM_KIND)                                 :: Tau2
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
!!$    IF( c%frill%leaf )THEN ! Leaf condition ? 
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
!!$    ENDIF
!!$
!!$    CALL SpAMM_redecorate_tree_2d_symm(c)
!!$
!!$  END SUBROUTINE SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur



end module spamm_nbdyalgbra_times
