module SpAMMsand_inverse_squareroot

  USE spammpack 
!  USE spammsand_structures
!  USE spammsand_rqi_extremals

  implicit none

contains

  SUBROUTINE spammsand_scaled_newton_shulz_inverse_squareroot( x, z, tau )

    TYPE(SpAMM_type_2d_symm) , POINTER, INTENT(IN)    :: x
    TYPE(SpAMM_type_2d_symm) , POINTER, INTENT(INOUT) :: z
    REAL(SpAMM_KIND),                   INTENT(IN)    :: Tau

    TYPE(SpAMM_type_2d_symm) , POINTER                :: x =>NULL(),t=>NULL()
    INTEGER                                           :: i,j
    REAL(SpAMM_KIND)                                  :: EvMin, EvMax, sc, xo, xo_prev 
    REAL(SpAMM_KIND)                                  :: xo_analytic, xmax,delta, FillN, FillN_prev

    s%frill%ndimn
    !
    z => SpAMM_identity_matrix(M,N)
    x => SpAMM_identity_matrix(M,N)
    t => SpAMM_identity_matrix(M,N)
    !
    FillN=1d10
    !
    DO i = 1, 22

       ! |X_n> = <Z_n|S> |Z_n>
       t => SpAMM_tree_2d_symm_t_times_tree_2d_symm( z, s, tau*1d-3, in=t1, alpha=SpAMM_zero, beta=SpAMM_one)
       x => SpAMM_tree_2d_symm_n_times_tree_2d_symm( t, z, tau     , in=x,  alpha=SpAMM_zero, beta=SpAMM_one)
       
       ! monitor the trace for convergence, maybe look at rate of change at some point too?:
       FillN_prev=FillN
       FillN=DBLE(N-Trace(X))/DBLE(N)       
       !        
       IF(FillN>0.4d0)THEN
          delta=1d-1  ! maybe this should be a variable too, passed in?
          X => GSOLVE_SPECTRAL_Shift( X, low_prev=0d0, high_prev=1d0, low_new=delta, high_new=1d0-delta )
          sc=ScaleInvSqrt(0d0)
       ELSE
          sc=1d0
       ENDIF

       X =>  GSOLVE_SPECTRAL_InvSqrt_Scaled_NS( X, sc ) 

       ! |Z_n+1> =  <Z_n| X_n>  n or t here ?
!       t => SpAMM_tree_2d_symm_n_times_tree_2d_symm( z, x, tau, in=t, alpha=SpAMM_zero, beta=SpAMM_one)
       t => SpAMM_tree_2d_symm_t_times_tree_2d_symm( z, x, tau, in=t, alpha=SpAMM_zero, beta=SpAMM_one)
       ! unfortunately, the update could not be done in place 
       CALL SpAMM_tree_2d_symm_copy_tree_2d_symm(z,t)

       ! end debug
       !--------------------------------------------------------------------------------

       IF(FillN<0.1d0.AND.FillN>FillN_prev)THEN
          RETURN
       ELSEIF(FillN<Tau)THEN
          RETURN
`       ENDIF

       ! best acceleration we can hope for
       xo_analytic=xo_analytic*(9d0/4d0)*sc
       !

    ENDDO
   
    call delete(x)
    call delete(t)

  END SUBROUTINE  GSOLVE_SCALED_NEWTON_SCHULZ_INVSQT
