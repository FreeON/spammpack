module spamm_nbdyalgbra_times 

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration
  use spamm_elementals

  implicit none

CONTAINS
  !++NBODYTIMES: SpAMM generalized n-body algebras for times ____ NBODYTIMES _________________
  !++NBODYTIMES: generalized products, dots, contractions & convolutions (X)
  !++NBODYTIMES:   ... [TREE-ONE-D X TREE-ONE-D] ... [TREE-ONE-D X TREE-ONE-D] ...   
  !++NBODYTIMES:   SpAMM_tree_1d_dot_tree_1d_recur
  !++NBODYTIMES:     dot = (a_1,b_1) 
  RECURSIVE FUNCTION SpAMM_tree_1d_dot_tree_1d_recur(a, b ) result(dot)

    TYPE(SpAMM_tree_1d), POINTER :: a,b
    REAL(SpAMM_KIND)             :: Dot

    Dot=SpAMM_Zero

    if(.not.associated(a))return
    if(.not.associated(b))return

    if(a%frill%leaf)then
       dot=DOT_PRODUCT(a%chunk(1:SBS),b%chunk(1:SBS))
    else
       dot = dot + SpAMM_tree_1d_dot_tree_1d_recur(a%child_0, b%child_0)
       dot = dot + SpAMM_tree_1d_dot_tree_1d_recur(a%child_1, b%child_1)
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

    !    CALL SpAMM_print_tree_1d_recur (randm) 

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
  FUNCTION SpAMM_tree_2d_symm_times_tree_1d(a, b, Tau, alpha_O, beta_O, c) RESULT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A
    TYPE(SpAMM_tree_1d),      POINTER,           INTENT(IN)    :: B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
    TYPE(SpAMM_tree_1d),      POINTER, OPTIONAL, INTENT(INOUT) :: C
    TYPE(SpAMM_tree_1d),      POINTER                          :: D
    REAL(SpAMM_KIND)                                           :: alpha, beta
    INTEGER                                                    :: Depth
    REAL(SpAMM_KIND)                                           :: Tau2

    ! figure the starting conditions ...
    if(present(c))then
       d => c
    else
       d => NULL()
    endif
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    ! need a new tree top? lets instantiate one ... 
    if(.not.associated(d))&
       d => SpAMM_new_top_tree_1d(a%frill%ndimn(1))

    if(present(alpha_O))then
       d=>SpAMM_scalar_times_tree_1d(alpha_O, d)
    endif
    
    beta =SpAMM_one
    if(present( beta_O))beta = beta_O

    Depth=0
    CALL SpAMM_tree_2d_symm_times_tree_1d_recur(d, A, beta, B, Tau2, Depth )

  END FUNCTION SpAMM_tree_2d_symm_times_tree_1d

  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_1d_recur(C, A, beta, B, Tau2, Depth ) !<++NBODYTIMES|
   !                    c_1 => alpha*c_1 + beta*(aT_2.b_1) (recursive)                     !<++NBODYTIMES|    

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A
    TYPE(SpAMM_tree_1d),      POINTER, INTENT(IN)    :: B
    TYPE(SpAMM_tree_1d),      POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_1d), POINTER                     :: c0,c1

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return 

    IF( c%frill%leaf )THEN ! Leaf condition ? 

!       write(*,33) c%frill%bndbx, a%frill%bndbx, b%frill%bndbx
!33     format("[",i2,"-",i2,"]=[",i2,"-",i2,"]x[",i2,"-",i2,"] . [",i2,"=",i2,"]")

       c%chunk(1:SBS)=c%chunk(1:SBS)+beta*MATMUL(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS))
       c%frill%flops  = c%frill%flops + SBS2 + 2*SBS

!       WRITE(*,44)c%chunk
44     format(16(F10.6,", "))

    ELSE

       ! children's place on the tree-1d
       c0=> SpAMM_construct_tree_1d_0( c )
       c1=> SpAMM_construct_tree_1d_1( c )

       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c0, a%child_00, beta, b%child_0, Tau2, Depth+1 )      
       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c1, a%child_10, beta, b%child_0, Tau2, Depth+1 )      
       
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c0, a%child_01, beta, b%child_1, Tau2,  Depth+1 )             
       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c1, a%child_11, beta, b%child_1, Tau2,  Depth+1 )      
       !
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

    depth=0
    d=>a
    if(.not.associated(a))return

    CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, d)

  END FUNCTION SpAMM_scalar_times_tree_2d_symm

  !++NBODYTIMES:     SpAMM_scalar_times_tree_2d_recur
  !++NBODYTIMES:       a_1 => alpha*a_1 (recursive)
  recursive subroutine SpAMM_scalar_times_tree_2d_symm_recur(alpha, a)

    type(SpAMM_tree_2d_symm), pointer :: a
    real(SpAMM_KIND)                  :: alpha

    if(.not.associated(a))return

    IF(a%frill%leaf)THEN

       a%chunk(1:SBS,1:SBS)=alpha*a%chunk(1:SBS,1:SBS)
       a%frill%flops=a%frill%flops+SBS

    ELSE

       ! child along [00]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_00 ) !,depth+1 )
       ! child along [01]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_01 ) !,depth+1 )
       ! child along [10]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_10 ) !,depth+1 )
       ! child along [11]: 
       CALL SpAMM_scalar_times_tree_2d_symm_recur(alpha, a%child_11 ) !,depth+1 )

    ENDIF  

    ! merge & regarnish back up the tree ...
    CALL SpAMM_redecorate_tree_2d_symm(a)
    !
  end subroutine SpAMM_scalar_times_tree_2d_symm_recur


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (wrapper)
  FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau, NT_O, alpha_O, beta_O, in_O) RESULT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    LOGICAL, OPTIONAL,                           INTENT(IN)    :: NT_O
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: In_O
    TYPE(SpAMM_tree_2d_symm), POINTER                          :: D
    REAL(SpAMM_KIND)                                           :: alpha, beta
    INTEGER                                                    :: Depth
    LOGICAL                                                    :: NT
    REAL(SpAMM_KIND)                                           :: Tau2

    ! figure the starting conditions ...
    if(present(in_O))then
       d => in_O
    else
       d => NULL()
    endif
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    ! need a new tree? then instantiate one ... 
    if(.not.associated(d))&
       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn)

    if(present(NT_O))then
       NT=NT_O    ! If NT_O==FALSE, then A^t.B
    else
       NT=.TRUE.  ! default is A.B
    endif

    if(present(alpha_O))then
       d => SpAMM_scalar_times_tree_2d_symm(alpha_O, d)
    endif

    beta =SpAMM_one
    if(present( beta_O))beta = beta_O

    CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(d, A, NT, B, Tau2, Depth, beta )

  END FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur(C, A, NT, B, Tau2, Depth,  beta )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    LOGICAL,                           INTENT(IN)    :: NT
    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c10,c11

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return  

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       IF(NT)THEN
          c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+beta*MATMUL(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))
       ELSE
          c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+ &
               beta*MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
       ENDIF
       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2

   ELSE

       ! find some memory ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c10=>SpAMM_construct_tree_2d_symm_10(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first pass ...
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c00,a%child_00,NT,b%child_00,Tau2,Depth+1,beta)     
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c01,a%child_00,NT,b%child_01,Tau2,Depth+1,beta)      

       IF(NT)THEN
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c10,a%child_10,NT,b%child_00,Tau2,Depth+1,beta)      
       ELSE
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c10,a%child_01,NT,b%child_00,Tau2,Depth+1,beta)      
       ENDIF

       IF(NT)THEN
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c11,a%child_10,NT,b%child_01,Tau2,Depth+1,beta)      
       ELSE
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c11,a%child_01,NT,b%child_01,Tau2,Depth+1,beta)      
       ENDIF

       ! ... & another pass 
       IF(NT)THEN
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c00,a%child_01,NT,b%child_10,Tau2,Depth+1,beta)      
       ELSE
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c00,a%child_10,NT,b%child_10,Tau2,Depth+1,beta)      
       ENDIF

       IF(NT)THEN
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c01,a%child_01,NT,b%child_11,Tau2,Depth+1,beta)      
       ELSE
          CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c01,a%child_10,NT,b%child_11,Tau2,Depth+1,beta)      
       ENDIF

       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c10,a%child_11,NT,b%child_10,Tau2,Depth+1,beta)      
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c11,a%child_11,NT,b%child_11,Tau2,Depth+1,beta)      
       !
    ENDIF

    CALL SpAMM_redecorate_tree_2d_symm(c)

  END SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur












  FUNCTIOn SpAMM_tree_2d_symm_T_times_tree_2d_symm(a, b, Tau, alpha_O, beta_O, in_O) RESULT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: In_O
    TYPE(SpAMM_tree_2d_symm), POINTER                          :: D
    REAL(SpAMM_KIND)                                           :: alpha, beta
    INTEGER                                                    :: Depth
    REAL(SpAMM_KIND)                                           :: Tau2
    ! figure the starting conditions ...
    if(present(in_O))then
       d => in_O
    else
       d => NULL()
    endif
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    ! need a new tree? then instantiate one ... 
    if(.not.associated(d))&
       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn)

    if(present(alpha_O))then
       d => SpAMM_scalar_times_tree_2d_symm(alpha_O, d)
    endif

    beta =SpAMM_one
    if(present( beta_O))beta = beta_O

    CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(d, A, B, Tau2, Depth, beta )

  END FUNCTION SpAMM_tree_2d_symm_T_times_tree_2d_symm


  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(C, A, B, Tau2, Depth,  beta )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c10,c11

    integer, dimension(1:2) :: ahi,alo , bhi,blo , chi , clo, athi , atlo

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return  

    IF( c%frill%leaf )THEN ! Leaf condition ? 
       
!       WRITE(*,*)'before ',SQRT(SUM(c%chunk**2)),' beta = ',beta,SQRT(SUM(a%chunk**2)),SQRT(SUM(b%chunk**2))

 !      blo=b%frill%bndbx(0,:)
 !      bhi=b%frill%bndbx(1,:)
 !      clo=c%frill%bndbx(0,:)
 !      chi=c%frill%bndbx(1,:)
 !      alo=a%frill%bndbx(0,:)
 !      ahi=a%frill%bndbx(1,:)

 !       aTlo(1)=a%frill%bndbx(0,2)
 !      aTlo(2)=a%frill%bndbx(0,1)
 !      aThi(1)=a%frill%bndbx(1,2)
 !      aThi(2)=a%frill%bndbx(1,1)


!       WRITE(*,67)clo(1),chi(1), clo(2),chi(2), &  
!                  aTlo(1),aThi(1), aTlo(2),aThi(2), & 
!                  blo(1),bhi(1), blo(2),bhi(2)

!       WRITE(*,77)chi(1)-clo(1)+1, chi(2)-clo(2)+1, &  
!                  ahi(1)-alo(1)+1, ahi(2)-alo(2)+1, & 
!                  bhi(1)-blo(1)+1, bhi(2)-blo(2)+1

67     format('[',I3,'-',I3,', ',I3,'-',I3,'] = [',I3,'-',I3,', ',I3,'-',I3,']^t x [',I3,'-',I3,', ',I3,'-',I3,']')
77     format('[',I3,', ',I3,'] = [',I3,', ',I3,'] x [',I3,', ',I3,']')

!       c%chunk(1:(chi(1)-clo(1)+1), 1:(chi(2)-clo(2)+1)) = &
!       c%chunk(1:(chi(1)-clo(1)+1), 1:(chi(2)-clo(2)+1)) + &
!       beta*MATMUL( TRANSPOSE( a%chunk( 1:ahi(1)-alo(1)+1 , 1:ahi(2)-alo(2)+1 ) ), &
!                   b%chunk(1:(bhi(1)-blo(1)+1), 1:(bhi(2)-blo(2)+1)))

       c%chunk(1:SBS,1:SBS)=c%chunk(1:SBS,1:SBS)+ &
               beta*MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))

!       WRITE(*,*)'after ',SQRT(SUM(c%chunk**2))

       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2

   ELSE

       ! find some memory ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c10=>SpAMM_construct_tree_2d_symm_10(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first pass ...
!       WRITE(*,*)' 00 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c00, a%child_00, b%child_00, & 
                                             Tau2, Depth+1, beta )     

!       WRITE(*,*)' 01=00:01 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c01, a%child_00, b%child_01, & 
                                             Tau2, Depth+1, beta )      

!       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c10, a%child_10, b%child_00, & 
!                                             Tau2, Depth+1, beta )      

 !      WRITE(*,*)' 10=01^t:00 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c10, a%child_01, b%child_00, & 
                                             Tau2, Depth+1, beta )      

!       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c11, a%child_10, b%child_01, & 
!                                             Tau2, Depth+1, beta )      

!       WRITE(*,*)' 11=01^t:01 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c11, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, beta )      

!       WRITE(*,*)' ------------------'

       ! ... & another pass 
!       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c00, a%child_01, b%child_10, & 
!                                             Tau2, Depth+1, beta )      
!       WRITE(*,*)' 00=10^t:10 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c00, a%child_10, b%child_10, & 
                                             Tau2, Depth+1, beta )      


!       WRITE(*,*)' 01=10^t:11 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c01, a%child_10, b%child_11, & 
                                             Tau2, Depth+1, beta )      

!       WRITE(*,*)' 10=11:10 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c10, a%child_11, b%child_10, & 
                                             Tau2, Depth+1, beta )      

!       WRITE(*,*)' 11=11:11 '
       CALL SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur(c11, a%child_11, b%child_11, & 
                                             Tau2, Depth+1, beta )      

       !
    ENDIF

    CALL SpAMM_redecorate_tree_2d_symm(c)

  END SUBROUTINE SpAMM_tree_2d_symm_T_times_tree_2d_symm_n_recur



end module spamm_nbdyalgbra_times
