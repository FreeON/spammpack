module spamm_nbdyalgbra_times 

  use spamm_structures
  use spamm_xstructors
  use spamm_decoration

  implicit none

CONTAINS
  !++NBODYTIMES: SpAMM generalized n-body algebras for times ____ NBODYTIMES _________________
  !++NBODYTIMES: generalized products, dots, contractions & convolutions (X)
  !++NBODYTIMES:   ... [TREE-ONE-D X TREE-ONE-D] ... [TREE-ONE-D X TREE-ONE-D] ...   
  !++NBODYTIMES:   SpAMM_tree_1d_dot_tree_1d_recur
  !++NBODYTIMES:     dot = (a_1,b_1) 
  RECURSIVE FUNCTION SpAMM_tree_1d_dot_tree_1d_recur(a, b ) result(dot)

    TYPE(SpAMM_tree_1d), POINTER :: a,b
    REAL(SpAMM_KIND)             :: Dot, Dot0, Dot1

    Dot=SpAMM_Zero

    if(.not.associated(a))return
    if(.not.associated(b))return

    if(a%leaf)then
       dot=DOT_PRODUCT(a%chunk(1:SBS),b%chunk(1:SBS))
    else
       dot = dot + SpAMM_tree_1d_dot_tree_1d_recur(a%child_0, b%child_0)
       dot = dot + SpAMM_tree_1d_dot_tree_1d_recur(a%child_1, b%child_1)
    endif
  END FUNCTION SpAMM_tree_1d_dot_tree_1d_recur
  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-ONE-D] ... [TREE-TWO-D X TREE-ONE-D] ...   
  !++NBODYTIMES:     SpAMM_tree_2d_symm_times_tree_1d
  !++NBODYTIMES:     c_1 => alpha*c_1 + beta*(a_2.b_1) (wrapper)
  FUNCTION SpAMM_tree_2d_symm_times_tree_1d(a, b, Tau, alpha, beta, c) RESUT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A
    TYPE(SpAMM_tree_1d),      POINTER,           INTENT(IN)    :: B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
    TYPE(SpAMM_tree_1d,       POINTER, OPTIONAL, INTENT(INOUT) :: C
    TYPE(SpAMM_tree_1d,       POINTER                          :: D
    REAL(SpAMM_KIND)                                           :: alpha, beta
    INTEGER                                                    :: Depth
    REAL(SpAMM_KIND)                                           :: Tau2

    ! figure the starting conditions ...
    if(present(c))then
       d => c
    else
       d => NULL()
       if(present(alpha))write(*,*) ' multiplying though by ',alpha, ' whilst C un-inited ...'
    endif
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    ! need a new tree (not there, or a  temp), lets instantiate one ... 
    if(.not.associated(d))&
       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn(1))

    alpha=SpAMM_one
    beta =SpAMM_one
    if(present(alpha_O))alpha=alpha_O
    if(present( beta_O))beta = beta_O

    Depth=0
    CALL SpAMM_tree_2d_symm_times_tree_1d_recur(d, A, B, Tau2, Depth, alpha, beta )

  END FUNCTION SpAMM_tree_2d_symm_times_tree_1d

  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_1d_symm_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (recursive)










  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_1d_recur(C, A, B, Tau2, Depth, alpha, beta )

    TYPE(SpAMM_tree_1d), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_1d), TARGET,  INTENT(INOUT) :: C
    REAL(SpAMM_KIND),             INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                            :: Tau2
    INTEGER                                     :: Depth
    TYPE(SpAMM_tree_1d), POINTER                :: c0,c1

    if(.not.associated(a))return
    if(.not.associated(b))return

    if(a%Norm2*b%Norm2<=Tau2)return    ! n-body for matrices with decay (and some structure)

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS) = alpha*c%chunk(1:SBS) + beta*matmul(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS))
       c%frill%flops  = c%frill%flops + SBS2 + 2*SBS

    ELSE

       ! find some memory ...
       c0=>SpAMM_construct_tree_1d_symm_0(c)
       c1=>SpAMM_construct_tree_1d_symm_1(c)

       ! a first pass with [0], [1] memory ...
       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c0, a_child_00, b%child_0, & 
                                      Tau2, Depth+1, alpha=alpha, beta=beta )      

       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c1, a_child_10, b%child_0, & 
                                      Tau2, Depth+1, alpha=alpha, beta=beta )      

       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c0, a_child_01, b%child_1, & 
                                      Tau2, Depth+1, alpha=alpha, beta=beta )      

       CALL SpAMM_tree_2d_symm_times_tree_1d_recur(c1, a_child_11, b%child_1, & 
                                      Tau2, Depth+1, alpha=alpha, beta=beta )      
       !
    ENDIF
  END SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_recur

  !++NBODYTIMES:   ... [TREE-TWO-D X TREE-TWO-D] ... [TREE-TWO-D X TREE-TWO-D] ...   
  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (wrapper)
  FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm(a, b, Tau, alpha, beta, c) RESUT(d)

    TYPE(SpAMM_tree_2d_symm), POINTER,           INTENT(IN)    :: A, B
    REAL(SpAMM_KIND),                            INTENT(IN)    :: Tau
    REAL(SpAMM_KIND),                  OPTIONAL, INTENT(IN)    :: alpha_O, beta_O
    TYPE(SpAMM_tree_1d,       POINTER, OPTIONAL, INTENT(INOUT) :: C
    TYPE(SpAMM_tree_1d,       POINTER                          :: D
    REAL(SpAMM_KIND)                                           :: alpha, beta
    INTEGER                                                    :: Depth
    REAL(SpAMM_KIND)                                           :: Tau2

    ! figure the starting conditions ...
    if(present(c))then
       d => c
    else
       d => NULL()
       if(present(alpha))write(*,*) ' multiplying though by ',alpha, ' whilst C un-inited ...'
    endif
    ! bail if we can ...
    if(.not.associated(a))return
    if(.not.associated(b))return

    ! here is the squared threshold 
    Tau2=Tau*Tau

    ! need a new tree (not there, or a  temp), lets instantiate one ... 
    if(.not.associated(d))&
       d => SpAMM_new_top_tree_2d_symm(a%frill%ndimn(1))

    alpha=SpAMM_one
    beta =SpAMM_one
    if(present(alpha_O))alpha=alpha_O
    if(present( beta_O))beta = beta_O

    Depth=0
    CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(d, A, B, Tau2, Depth, alpha, beta )

  END FUNCTION SpAMM_tree_2d_symm_times_tree_2d_symm

  !++NBODYTIMES:   SpAMM_tree_2d_symm_times_tree_2d_symm_recur
  !++NBODYTIMES:     c_2 => alpha*c_2 + beta*(a_2.b_2) (recursive)
  RECURSIVE SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur(C, A, B, Tau2, Depth, alpha, beta )

    TYPE(SpAMM_tree_2d_symm), POINTER, INTENT(IN)    :: A, B
    TYPE(SpAMM_tree_2d_symm), TARGET,  INTENT(INOUT) :: C
    REAL(SpAMM_KIND),                  INTENT(IN)    :: alpha, beta
    REAL(SpAMM_KIND)                                 :: Tau2
    INTEGER                                          :: Depth
    TYPE(SpAMM_tree_2d_symm), POINTER                :: c0,c1

    if(.not.associated(a))return
    if(.not.associated(b))return

    if(a%Norm2*b%Norm2<=Tau2)return    ! n-body for matrices with decay (and some structure)

    IF( c%frill%leaf )THEN ! Leaf condition ? 

       c%chunk(1:SBS:1:SBS)=alpha*c%chunk(1:SBS:SBS)+beta*matmul(a%chunk(1:SBS,1:SBS),b%chunk(1:SBS,1:SBS))
       c%frill%flops = c%frill%flops + SBS3 + 2*SBS2

    ELSE

       ! find some memory ...
       c00=>SpAMM_construct_tree_2d_symm_00(c)
       c01=>SpAMM_construct_tree_2d_symm_01(c)
       c11=>SpAMM_construct_tree_2d_symm_11(c)

       ! a first pass with [00], [01] & [11] memory ...
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c00, a_child_00, b%child_00, & 
                                             Tau2, Depth+1, alpha=alpha, beta=beta )      
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c01, a_child_00, b%child_01, & 
                                             Tau2, Depth+1, alpha=alpha, beta=beta )      
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c11, a_child_10, b%child_01, & 
                                             Tau2, Depth+1, alpha=alpha, beta=beta )      
       ! ... & another pass 
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c00, a_child_01, b%child_10, & 
                                             Tau2, Depth+1, alpha=alpha, beta=beta )      
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c01, a_child_01, b%child_11, & 
                                             Tau2, Depth+1, alpha=alpha, beta=beta )      
       CALL SpAMM_tree_2d_symm_times_tree_2d_symm_recur(c11, a_child_11, b%child_11, & 
                                             Tau2, Depth+1, alpha=alpha, beta=beta )      
       !
    ENDIF

  END SUBROUTINE SpAMM_tree_2d_symm_times_tree_2d_symm_recur

end module spamm_nbdyalgbra_times
