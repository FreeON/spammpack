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

    WRITE(*,*)'A randm = ',randm%frill%bndbx

    depth=0
!    CALL init_random_seed()    

    CALL SpAMM_random_unormalized_tree_1d_recur (randm, depth)

    WRITE(*,*)'B randm = ',randm%frill%bndbx
    STOP

    CALL SpAMM_print_tree_1d_recur (randm) 

    ! normalize the vector ...
    renorm=SpAMM_one/sqrt(randm%frill%norm2)

    WRITE(*,*)' b '

    WRITE(*,*)' RENORM =',RENORM,' randm = ',randm%frill%bndbx
    randm=>SpAMM_scalar_times_tree_1d(renorm, randm)

    WRITE(*,*)' c '


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

    if(.not.associated(randm))return

    IF(randm%frill%leaf)THEN
       CALL RANDOM_NUMBER(randm%chunk)
       write(*,*)randm%frill%bndbx,' chunk = ',randm%chunk
    ELSE

       WRITE(*,*)'A depth = ',depth,' rand%bb = ',randm%frill%bndbx
       ! child along [0]: [lo,mid] ... 
       CALL SpAMM_random_unormalized_tree_1d_recur(SpAMM_construct_tree_1d_0(randm), depth+1 )
       WRITE(*,*)'B depth = ',depth,' rand%bb = ',randm%frill%bndbx
       ! child along [1]: [mid+1,hi] ... 
       CALL SpAMM_random_unormalized_tree_1d_recur(SpAMM_construct_tree_1d_1(randm), depth+1 )
       WRITE(*,*)'C depth = ',depth,' rand%bb = ',randm%frill%bndbx
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

    CALL SpAMM_scalar_times_tree_1d_recur(alpha, d,depth)

  END FUNCTION SpAMM_scalar_times_tree_1d

  !++NBODYTIMES:     SpAMM_scalar_times_tree_1d_recur
  !++NBODYTIMES:       a_1 => alpha*a_1 (recursive)


  recursive subroutine SpAMM_scalar_times_tree_1d_recur(alpha, a, depth)

    type(SpAMM_tree_1d), pointer :: a
    real(SpAMM_KIND)             :: alpha

integer :: depth

    if(.not.associated(a))return

    write(*,*)depth,a%frill%leaf,a%frill%bndbx


    if(depth>6)stop

    IF(a%frill%leaf)THEN

       a%chunk=alpha*a%chunk
       a%frill%flops=a%frill%flops+SBS

    ELSE
       ! child along [0]: [lo,mid] ... 
       CALL SpAMM_scalar_times_tree_1d_recur(alpha, SpAMM_construct_tree_1d_0(a),depth+1 )
       ! child along [1]: [mid+1,hi] ... 
       CALL SpAMM_scalar_times_tree_1d_recur(alpha, SpAMM_construct_tree_1d_1(a),depth+1 )
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

    alpha=SpAMM_one
    beta =SpAMM_one
    if(present(alpha_O))alpha=alpha_O
    if(present( beta_O))beta = beta_O

    Depth=0
    CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(d, A, B, Tau2, Depth, alpha, beta )

  END FUNCTION SpAMM_tree_2d_symm_times_tree_1d

  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_n_times_tree_1d_recur(C, A, B, Tau2, Depth, alpha, beta ) !<++NBODYTIMES|
  !                    c_1 => alpha*c_1 + beta*(a_2.b_1) (recursive)                              !<++NBODYTIMES|    

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A
    TYPE(SpAMM_tree_1d),      POINTER, INTENT(IN)    :: B
    TYPE(SpAMM_tree_1d),      POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_1d), POINTER                     :: c0,c1

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return 

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS) = alpha*c%chunk(1:SBS) + beta*matmul(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS))
       c%frill%flops  = c%frill%flops + SBS2 + 2*SBS

    ELSE

       ! find some memory ...
       c0=>SpAMM_construct_tree_1d_0(c)
       c1=>SpAMM_construct_tree_1d_1(c)

       ! a first pass with [0], [1] memory ...
       CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(c0, a%child_00, b%child_0, & 
                                      Tau2, Depth+1, alpha, beta )      

       CALL SpAMM_tree_2d_symm_t_times_tree_1d_recur(c1,a%child_01, b%child_0, & 
                                      Tau2, Depth+1, alpha, beta )      
       
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(c0, a%child_01, b%child_1, & 
            Tau2, Depth+1, alpha, beta )      
       
       CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(c1, a%child_11, b%child_1, & 
                                      Tau2, Depth+1, alpha, beta )      
       !
    ENDIF
  END SUBROUTINE SpAMM_tree_2d_symm_n_times_tree_1d_recur


  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_t_times_tree_1d_recur(C, A, B, Tau2, Depth, alpha, beta ) !<++NBODYTIMES|
   !                    c_1 => alpha*c_1 + beta*(aT_2.b_1) (recursive)                              !<++NBODYTIMES|    

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A
    TYPE(SpAMM_tree_1d),      POINTER, INTENT(IN)    :: B
    TYPE(SpAMM_tree_1d),      POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_1d), POINTER                     :: c0,c1

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return 

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS) = alpha*c%chunk(1:SBS) + beta*MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS))
       c%frill%flops  = c%frill%flops + SBS2 + 2*SBS

    ELSE

       ! find some memory ...
       c0=>SpAMM_construct_tree_1d_0(c)
       c1=>SpAMM_construct_tree_1d_1(c)

       ! a first pass with [0], [1] memory ...
       CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(c0, a%child_00, b%child_0, & 
                                      Tau2, Depth+1, alpha, beta )      

       CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(c1,a%child_01, b%child_0, & 
                                      Tau2, Depth+1, alpha, beta )      
       
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_t_times_tree_1d_recur(c0, a%child_01, b%child_1, & 
            Tau2, Depth+1, alpha, beta )      
       
       CALL SpAMM_tree_2d_symm_n_times_tree_1d_recur(c1, a%child_11, b%child_1, & 
                                      Tau2, Depth+1, alpha, beta )      
       !
    ENDIF
  END SUBROUTINE SpAMM_tree_2d_symm_t_times_tree_1d_recur



  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-TWO-D] ... [TREE-TWO-D X TREE-TWO-D] ...   
  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (wrapper)
  FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau, alpha_O, beta_O, c_O, LeftT_O, RghtT_O) RESULT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
    TYPE(SpAMM_tree_2d_symm), POINTER, OPTIONAL, INTENT(INOUT) :: C_O
    LOGICAL,                           OPTIONAL                :: LeftT_O,RghtT_O
    TYPE(SpAMM_tree_2d_symm), POINTER                          :: D
    REAL(SpAMM_KIND)                                           :: alpha, beta
    INTEGER                                                    :: Depth
    REAL(SpAMM_KIND)                                           :: Tau2
    LOGICAL                                                    :: LeftT,RghtT
    ! figure the starting conditions ...
    if(present(c_O))then
       d => c_O
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

    alpha=SpAMM_one
    beta =SpAMM_one
    if(present(alpha_O))alpha=alpha_O
    if(present( beta_O))beta = beta_O

    LeftT=.FALSE.
    RghtT=.FALSE.
    if(present(LeftT_O))LeftT=LeftT_O
    if(present(RghtT_O))RghtT=RghtT_O
    IF(LeftT.AND.RghtT)STOP ' dumbass '

    Depth=0
    IF(LeftT)THEN
       CALL SpAMM_tree_2d_symm_t_times_tree_2d_symm_n_recur(d, A, B, Tau2, Depth, alpha, beta )
    ELSEIF(RghtT)THEN
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(d, A, B, Tau2, Depth, alpha, beta )
    ELSE
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(d, A, B, Tau2, Depth, alpha, beta )
    ENDIF


  END FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm

  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(C, A, B, Tau2, Depth, alpha, beta )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c11

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return  

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS,1:SBS)=alpha*c%chunk(1:SBS,1:SBS) + &
            beta*MATMUL(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))

       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2

    ELSE

       ! find some memory ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first pass with [00], [01] & [11] memory ...
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c00, a%child_00, b%child_00, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c01, a%child_00, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_t_times_tree_2d_symm_n_recur(c11, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(c00, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c01, a%child_01, b%child_11, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c11, a%child_11, b%child_11, & 
                                             Tau2, Depth+1, alpha, beta )      
       !
    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur

  !++NBODYTIMES:   SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.bT_2) (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(C, A, B, Tau2, Depth, alpha, beta )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c11

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return  

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS,1:SBS)=alpha*c%chunk(1:SBS,1:SBS) +  &
            beta*MATMUL(a%chunk(1:SBS,1:SBS),TRANSPOSE(b%chunk(1:SBS,1:SBS)))
       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2

    ELSE

       ! find some memory ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first pass with [00], [01] & [11] memory ...
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c00, a%child_00, b%child_00, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(c01, a%child_00, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(c11, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(c00, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c01, a%child_01, b%child_11, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c11, a%child_11, b%child_11, & 
                                             Tau2, Depth+1, alpha, beta )      
       !
    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur


  !++NBODYTIMES:   SpAMM_tree_2d_symm_t_times_tree_2d_symm_n_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(aT_2.b_2) (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_t_times_tree_2d_symm_n_recur(C, A, B, Tau2, Depth, alpha, beta )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_2d_symm), POINTER                :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c00,c01,c11

    if(.not.associated(a))return
    if(.not.associated(b))return

    ! n-body occlusion & culling of the product for matrices with decay (and some structure)
    if(a%frill%Norm2*b%frill%Norm2<=Tau2)return  

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS,1:SBS)=alpha*c%chunk(1:SBS,1:SBS) + &
            beta*MATMUL(TRANSPOSE(a%chunk(1:SBS,1:SBS)),b%chunk(1:SBS,1:SBS))
       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2

    ELSE

       ! find some memory ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first pass with [00], [01] & [11] memory ...
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c00, a%child_00, b%child_00, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c01, a%child_00, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_t_times_tree_2d_symm_n_recur(c11, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_t_recur(c00, a%child_01, b%child_01, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c01, a%child_01, b%child_11, & 
                                             Tau2, Depth+1, alpha, beta )      
       CALL SpAMM_tree_2d_symm_n_times_tree_2d_symm_n_recur(c11, a%child_11, b%child_11, & 
                                             Tau2, Depth+1, alpha, beta )      

       !
    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_t_times_tree_2d_symm_n_recur



end module spamm_nbdyalgbra_times
