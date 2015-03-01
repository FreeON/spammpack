MODULE SpAMMSand_rqi_extremals

  USE spammpack 

CONTAINS
  !> SpAMM routines for spectral estimation (extremal eigenvalues)
  !!
  !! RQI for finding the min/max extrema of SpAMM matrix a:
  FUNCTION SpAMMSand_rqi_extremal(a,tau,high_O) RESULT(Omega)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: a
    REAL(SpAMM_KIND),      INTENT(IN)                 :: Tau
    LOGICAL, OPTIONAL                                 :: high_O

    TYPE(SpAMM_tree_1d), POINTER :: x=>NULL(),g=>NULL(),h=>NULL()
    TYPE(SpAMM_tree_1d), POINTER :: Ax=>NULL(),Ah=>NULL(),gOld=>NULL(),hOld=>NULL()
    INTEGER              :: CG, MinMax
    INTEGER, PARAMETER   :: NCG=100

    REAL(SpAMM_KIND)     :: omega,omega_old
    REAL(SpAMM_KIND)     :: xx, hh, xh, hx, xAx, xAh, hAx, hAh
    REAL(SpAMM_KIND)     :: sclr_Ax, sclr__x, beta, dot_g, dot_gold
    REAL(SpAMM_KIND)     :: LambdaPlus, LambdaMins, RQIPlus, RQIMins

    MINMAX=2 ! default is high (max) extremal
    IF(PRESENT(high_O))THEN
       IF(high_O)THEN
          MINMAX=2
       ELSE
          MINMAX=1
       ENDIF
    ENDIF

    IF(MinMax==2)THEN
       omega=-1d10      ! starting eigenvalue (to be maximized)
       gsign=-SpAMM_One ! sign of the cooresponding gradient
    ELSE
       omega= 1d10      ! starting eigenvalue (to be minimized)
       gsign=+SpAMM_One ! sign of the cooresponding gradient
    ENDIF

    M=a%frill%ndimn(1)              ! dimension of A along the [i] direction
    g   =>SpAMM_new_top_tree_1d(M)  ! gradient (analytic)
    h   =>SpAMM_new_top_tree_1d(M)  ! conjugate gradient (corrected g)
    Ax  =>SpAMM_new_top_tree_1d(M)  ! gradient with the matrix 
    Ah  =>SpAMM_new_top_tree_1d(M)  ! conjugate gradient with the matrix
    gOld=>SpAMM_new_top_tree_1d(M) 
    hOld=>SpAMM_new_top_tree_1d(M)

    x   =>SpAMM_random_tree_1d( M)  ! our extremal eigenvector
    CALL SpAMM_print_tree_1d_recur (x) 

    DO CG=1,NCG ! conjugate gradient iteration

       Ax => SpAMM_tree_2d_symm_times_tree_1d(a, x, Tau, alpha_O=SpAMM_zero, beta_O=SpAMM_one, c=Ax )
       xx =  SpAMM_tree_1d_dot_tree_1d_recur (x, x)
       xAx=  SpAMM_tree_1d_dot_tree_1d_recur (x,Ax)

       omega_old=omega
       omega=xAx/xx

       sclr_Ax = + gsign*SpAMM_Two/xx          ! g= + 2*(Ax-omega*x)/xx (minimizing)
       sclr__x = - gsign*SpAMM_Two*omega/xx    ! g= - 2*(Ax-omega*x)/xx (maximizing)

       g => SpAMM_tree_1d_plus_tree_1d ( g, sclr_Ax, Ax, inplace=SpAMM_zero) 
       g => SpAMM_tree_1d_plus_tree_1d ( g, sclr__x,  x)                     

       dot_g    = SpAMM_tree_1d_dot_tree_1d_recur( g,    g   )
       dot_gold = SpAMM_tree_1d_dot_tree_1d_recur( gOld, gOld)          

       IF(CG>1.AND.MOD(CG,15).NE.0)THEN
          IF(dot_gold/abs(omega).LE.1D-10)THEN
             ! if we are really close, steepest descents should be enuf ...
             beta=SpAMM_Zero
          ELSE
             beta=MAX(SpAMM_Zero,dot_g/dot_gold)
          ENDIF
       ELSE
          beta=SpAMM_Zero
       ENDIF

       ! convergence criteria
       IF( SQRT(dot_g)/ABS(Omega) < SQRT(Tau).AND.CG>16 &
            .OR. (MinMax==1.AND.Omega>Omega_old)        &
            .OR. (MinMax==2.AND.Omega<Omega_old) )EXIT

       ! the conjugate gradient, h = g + beta*hOld
       h => SpAMM_tree_1d_plus_tree_1d ( h, SpAMM_one, g, inplace=SpAMM_zero) 
       h => SpAMM_tree_1d_plus_tree_1d ( h,      beta, h)

       ! Ah = A.h

       Ah => SpAMM_tree_2d_symm_times_tree_1d(A, h, Tau, alpha_O=SpAMM_Zero, beta_O=SpAMM_One, c=Ah)

       hx =SpAMM_tree_1d_dot_tree_1d_recur(h,x)
       hh =SpAMM_tree_1d_dot_tree_1d_recur(h,h)
       xAh=SpAMM_tree_1d_dot_tree_1d_recur(x,Ah)
       hAx=SpAMM_tree_1d_dot_tree_1d_recur(h,Ax)
       hAh=SpAMM_tree_1d_dot_tree_1d_recur(h,Ah)        

       xh=hx      ! By symmetry
       hAx=xAh

       ! gOld = g; hOld = h
       gOld => SpAMM_tree_1d_copy_tree_1d (g, gOld) 
       hOld => SpAMM_tree_1d_copy_tree_1d (h, hOld) 

       ! roots of the line search (+/-) ...
       LambdaPlus=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx+SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
            -SpAMM_Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
            /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
       LambdaMins=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx-SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
            -SpAMM_Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
            /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))

       ! ... and cooresponding update of the objective
       RQIPlus=(xAx+LambdaPlus*(xAh+hAx)+hAh*LambdaPlus**2) &
            /( xx+LambdaPlus*(xh+hx)  +hh *LambdaPlus**2)
       RQIMins=(xAx+LambdaMins*(xAh+hAx)+hAh*LambdaMins**2) &
            /( xx+LambdaMins*(xh+hx)  +hh *LambdaMins**2)

       ! update of the eigenvector ..
       IF(MinMax==1)THEN ! for the minimizer ...          
          IF(RQIMins<RQIPlus)THEN
             x => SpAMM_tree_1d_plus_tree_1d(x, LambdaMins, h) ! x = x + LambdaMins*h
          ELSE             
             x => SpAMM_tree_1d_plus_tree_1d(x, LambdaPlus, h) ! x = x + LambdaPlus*h
          ENDIF
       ELSE ! for the maximizer ...
          IF(RQIMins>RQIPlus)THEN
             x => SpAMM_tree_1d_plus_tree_1d(x, LambdaMins, h) ! x = x + LambdaMins*h
          ELSE             
             x => SpAMM_tree_1d_plus_tree_1d(x, LambdaPlus, h) ! x = x + LambdaPlus*h
          ENDIF
       ENDIF

       IF(MinMax==1)THEN
          WRITE(*,33)omega,dot_g,CG
       ELSE
          WRITE(*,44)omega,dot_g,CG
       ENDIF

    END DO

33  FORMAT(' MIN E.V. = ',E24.16,', GRAD RQI = ',E16.8,' in ',I4,' NLCG steps')
44  FORMAT(' MAX E.V. = ',E24.16,', GRAD RQI = ',E16.8,' in ',I4,' NLCG steps')

    ! tidy ...
    CALL SpAMM_destruct_tree_1d_recur (x)
    CALL SpAMM_destruct_tree_1d_recur (g)
    CALL SpAMM_destruct_tree_1d_recur (h)
    CALL SpAMM_destruct_tree_1d_recur (Ax)
    CALL SpAMM_destruct_tree_1d_recur (Ah)
    CALL SpAMM_destruct_tree_1d_recur (gOld)
    CALL SpAMM_destruct_tree_1d_recur (hOld)

  END FUNCTION SpAMMSand_rqi_extremal

END MODULE SPAMMSAND_RQI_EXTREMALS 
