MODULE SpAMMSand_rqi_extremals
  USE SpAMM_nbdyalgebra
  USE SpAMM_conversion
  USE SpAMM_management
  USE SpAMM_types
  USE SpAMM_utilities
  USE test_utilities
CONTAINS
  !> SpAMM routines for spectral estimation (extremal eigenvalues)
  !!
  !! RQI for finding the min/max extrema of SpAMM matrix a:
  FUNCTION SpAMMSand_rqi_extremal(a,tau,high_O)

    TYPE(SpAMM_tree_2d_symm) , POINTER, INTENT(IN)    :: a
    REAL(SpAMM_KIND),      INTENT(IN)                 :: Tau

    TYPE(SpAMM_tree_1d), POINTER                      :: EVec=>NULL(),   x=>NULL(), &
                                                            g=>NULL(),   h=>NULL(), &
                                                           Ax=>NULL(),  Ah=>NULL(), &
                                                         xOld=>NULL(),gOld=>NULL(), &
                                                         hOld=>NULL()
    INTEGER              :: I,CG, MM
    INTEGER, PARAMETER   :: NCG=5000
  
    REAL(SpAMM_KIND)     :: beta,LambdaPlus,LambdaMins,RQIPlus,RQIMins,omega,omega_old
    REAL(SpAMM_KIND)     :: xx,hh,xh,hx,xAx,xAh,hAx,hAh,xnorm,dot_old

    MINMAX=2 ! default is high (max) extremal
    IF(PRESENT(high_O))THEN
       IF(high_O)THEN
          MINMAX=2
       ELSE
          MINMAX=1
       ENDIF
    ENDIF

    M=a%frill%ndimn(1) ! dimension of a along the [i] direction

    x   =>SpAMM_init_random_tree_1d(M) 
    g   =>SpAMM_new_top_tree_1d(M)
    h   =>SpAMM_new_top_tree_1d(M)
    Ax  =>SpAMM_new_top_tree_1d(M)
    Ah  =>SpAMM_new_top_tree_1d(M)
    xOld=>SpAMM_new_top_tree_1d(M)
    gOld=>SpAMM_new_top_tree_1d(M)
    hOld=>SpAMM_new_top_tree_1d(M)

    IF(MinMax==2)THEN
       omega=-1d10
       gsign= SpAMM_One
    ELSE
       omega= 1d10
       gsign=-SpAMM_One
    ENDIF

    DO CG=1,NCG

       Ax => SpAMM_tree_2d_symm_times_tree_1d(a, x)
       xx =  SpAMM_tree_1d_dot_tree_1d_recur (x, x)
       xAx=  SpAMM_tree_1d_dot_tree_1d_recur (x,Ax)

       omega_old=omega
       omega=xAx/xx

       sclr__x = + gsign*SpAMM_Two/xx
       sclr_Ax = - gsign*SpAMM_Two*omega/xx
       
       g => SpAMM_tree_1d_plus_tree_1d ( ax, x, scrl_Ax, sclr__x, g, SpAMM_Zero) 

       WRITE(*,*)'omega = ',omega,omega_old

       IF( SQRT(Dot(g,g)/ABS(Omega)) < CnvrgCrit .AND. CG>16 &
            .OR. (I==1.AND.Omega>Omega_old)                  &
            .OR. (I==2.AND.Omega<Omega_old)                    )EXIT
       
       IF(CG>1.AND.MOD(CG,15).NE.0)THEN
          dot_g    = SpAMM_tree_1d_dot_tree_1d_recur( g, g)
          dot_gold = SpAMM_tree_1d_dot_tree_1d_recur( gOld, gOld)
          IF(dot_gold.LE.1D-10)THEN
             beta=SpAMM_Zero
          ELSE
             beta=MAX(SpAMM_Zero,dot_g/dot_gold)
          ENDIF
        ELSE
          beta=SpAMM_Zero
        ENDIF

        ! h = g + beta*hOld
        h => SpAMM_tree_1d_plus_tree_1d ( g, hold, SpAMM_one, beta, h, SpAMM_zero) 
        ! Ah=A.h
        CALL Multiply(A%Root,h,Ah, Tau)

        hx =Dot(h,x)
        hh =Dot(h,h)
        xAh=Dot(x,Ah)
        hAx=Dot(h,Ax)
        hAh=Dot(h,Ah)

        ! By symmetry
        xh=hx
        hAx=xAh

        ! gOld = g; hOld = h
        gOld => SpAMM_tree_1d_copy_tree_1d (g, gOld) 
        hOld => SpAMM_tree_1d_copy_tree_1d (h, hOld) 

        LambdaPlus=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx+SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
          -SpAMM_Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
          /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
        LambdaMins=(SpAMM_Two*hh*xAx-SpAMM_Two*hAh*xx-SQRT((-SpAMM_Two*hh*xAx+SpAMM_Two*hAh*xx)**2     &
          -SpAMM_Four*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh)*(-(hx*xAx)-xAx*xh+SpAMM_Two*xAh*xx)))  &
          /(SpAMM_Two*(hAh*hx-SpAMM_Two*hh*xAh+hAh*xh))
        !
        RQIPlus=(xAx+LambdaPlus*(xAh+hAx)+hAh*LambdaPlus**2) &
          /( xx+LambdaPlus*(xh+hx)  +hh *LambdaPlus**2)
        RQIMins=(xAx+LambdaMins*(xAh+hAx)+hAh*LambdaMins**2) &
          /( xx+LambdaMins*(xh+hx)  +hh *LambdaMins**2)

        IF(I==1)THEN
          IF(RQIMins<RQIPlus)THEN
            ! x=x+LambdaMins*h
            CALL Add(x,SpAMM_One,h,LambdaMins)
          ELSE
            ! x=x+LambdaPlus*h
            CALL Add(x,SpAMM_One,h,LambdaPlus)
          ENDIF
        ELSE
          IF(RQIMins>RQIPlus)THEN
            ! x=x+LambdaMins*h
            CALL Add(x,SpAMM_One,h,LambdaMins)
          ELSE
            ! x=x+LambdaPlus*h
            CALL Add(x,SpAMM_One,h,LambdaPlus)
          ENDIF
        ENDIF

        IF(I==1)THEN
          WRITE(*,33)omega,SQRT(Dot(g,g))/ABS(Omega),CG
!          WRITE(*,33)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG
        ELSE
!          WRITE(*,44)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG
        ENDIF
      END DO

      IF(I==1)THEN
         WRITE(*,33)omega,SQRT(Dot(g,g))/ABS(Omega),CG
         WRITE(*,33)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG
         MinVec=>x
      ELSE
         WRITE(*,44)omega,SQRT(Dot(g,g))/ABS(Omega),CG
         WRITE(*,44)omega_dense,SQRT(Dot_product(g_dense,g_dense))/ABS(Omega_dense),CG         
         MaxVec=>x
      ENDIF

    ENDDO
!    WRITE(*,*)' EIGEN MIn = ',RQIMin,RQIMax
!    WRITE(*,*)' EIGEN MIn = ',RQIMin_dense,RQIMax_dense


33  FORMAT(' MIN E.V. = ',E24.16,', GRAD RQI = ',E16.8,' in ',I4,' NLCG steps')
44  FORMAT(' MAX E.V. = ',E24.16,', GRAD RQI = ',E16.8,' in ',I4,' NLCG steps')

    CALL Delete(g)
    CALL Delete(h)
    CALL Delete(Ax)
    CALL Delete(Ah)
    CALL Delete(xOld)
    CALL Delete(gOld)
    CALL Delete(hOld)

    CALL Delete(MinVec)
    CALL Delete(MaxVec)

    DEALLOCATE(x_dense)
    DEALLOCATE(g_dense)
    DEALLOCATE(h_dense)
    DEALLOCATE(Ax_dense)
    DEALLOCATE(Ah_dense)
    DEALLOCATE(xOld_dense)
    DEALLOCATE(gOld_dense)
    DEALLOCATE(hOld_dense)

  END SUBROUTINE GSOLVE_RQI_Extrema_RQI
END MODULE GSOLVE_RQI
