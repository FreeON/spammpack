module test_schulz

  REAL*8, PARAMETER ::  Approx3  = 2.85d00
  REAL*8, PARAMETER ::  ShiftSw  = 5.d-1
  REAL*8, PARAMETER ::  SPAMM_TAU= 1.d-11

contains

  REAL*8 function sq_scale(x,sc)
    REAL*8  :: x,sc 
    sq_scale=.5d0*SQRT(sc)*(3.0d0-sc*x)

!    WRITE(*,*)' x = ',x,' sq_scale = ',sq_scale

  end function sq_scale


  REAL*8 function rugged_sq(x,sc,delta)
    REAL*8  :: x,sc , delta
    rugged_sq=.5d0*SQRT(sc)*(3.0d0+delta-sc*(1d0-delta)*x)
  end function rugged_sq

  function minmax(x)
    REAL*8, DIMENSION(:) :: x
    REAL*8, DIMENSION(2) :: MinMax
    INTEGER :: I
    MinMax(1)=1d100
    MinMax(2)=-1d100
    DO I=1,SIZE(X)
       MinMax(1)=MIN(MinMax(1),X(I))
       MinMax(2)=MAX(MinMax(2),X(I))
    ENDDO

  end function minmax



  REAL*8 function sq(x)
    REAL*8  :: x
    sq=.5d0*(3.0d0-x)
  end function sq

  REAL*8 FUNCTION scale(xo)
    REAL*8 :: xo
    scale=MIN( Approx3, 3.d0/(1d0+SQRT(xo)+xo) )    
  END FUNCTION scale

  REAL*8 function shft(x, dx, low_old, high_old, low_new, high_new)
    REAL*8, INTENT(INOUT)  :: x, dx
    REAL*8, OPTIONAL :: low_old, high_old, low_new, high_new    
    dx=low_new-low_old
    shft=low_new+(x-low_old)*(high_new-low_new)/(high_old-low_old)
  END function shft




  SUBROUTINE SetScalarSpectrum(CondS,S,N)
    REAL*8 :: CondS
    REAL*8, DIMENSION(N) :: S
    integer,parameter :: seed = 86456
    REAL*8 :: low1,low2,high1,high2, dx, XMin,XMax
    INTEGER :: M

    goto 999


    call srand(seed)
    !    
    low1 = 1d10
    high1=-1d10
    DO i=1,N
       s(i)=RAND()
       low1 =MIN(low1, s(i))
       high1=MAX(high1,s(i))
    ENDDO         
    !
    low2=-LOG10(CondS)    
    high2=0

    DO i=1,N
       s(i)=shft(s(i), dx, low_old=low1, high_old=high1, low_new=low2, high_new=high2)
       s(i)=10d0**s(i)
!       WRITE(*,*)s(i)
    ENDDO


    RETURN

999 CONTINUE


    OPEN(UNIT=77,FILE='h2.evals')
!33_x8_11.evals')

    READ(77,*)M

    DO I=1,M
       READ(77,*)S(I)
    ENDDO

!    DO I=1,M
!       S(I)=S(I)/S(M)
!       WRITE(*,*)' SI = ',S(I)
!    ENDDO

    RETURN

    s=(/0.46231819888891973, &
        0.46231819888892006, &
        0.46947690964794725, & 
        0.46947690964794730, &
        0.83472877139968460, &
        0.83472877139968471, &
        0.86398915706764123, & 
        0.86398915706764134, &
        0.99985505910494166, &
        0.99985505910494199, &
        0.99985505910494221, & 
        0.99985505910494277, &
        1.0001449408950573,  &
        1.0001449408950573,  &
        1.0001449408950576,  & 
        1.0001449408950582,  & 
        1.1078193803996577,  & 
        1.1078193803996579,  &
        1.1379288929126554,  & 
        1.1379288929126559,  &
        1.5420058882896375,  &
        1.5420058882896377,  &
        1.5817328013938545,  & 
        1.5817328013938550 /)

    s=s/1.5817328013938550d0

  END SUBROUTINE SetScalarSpectrum

  SUBROUTINE ScaledNSHalf(N,S,Z,SpAMM_THRESHOLD)
    INTEGER              :: N
    REAL*8, DIMENSION(N) :: s,x,z,y,t,db1,db2
    REAL*8 :: SpAMM_THRESHOLD, dx, XMin,XMax, xo, xo_prev, sc, xo_analytic
    REAL*8 :: DIFF_ONE, DIFF_ISQ, DIFF_ONE_PREV, DIFF_ISQ_PREV,R1,R2
    LOGICAL :: NoScale

    DIFF_ONE=1d10
    DIFF_ISQ=1d10
    xo=1D10


    y=s
    WRITE(*,*)' S = ',MinMax(S)
    WRITE(*,*)' Z = ',MinMax(Z)
    WRITE(*,*) SpAMM_THRESHOLD

    NoScale=.FALSE.
    WRITE(*,31)
    DO j=1,30
       WRITE(*,*)' -------------------------',j,'----------------------------------'
       XMin= 1d10
       XMax=-1d10
       DO i=1,N
!          R1=2d0*(RAND()-0.5d0)*SpAMM_THRESHOLD
!          R2=2d0*(RAND()-0.5d0)*SpAMM_THRESHOLD

!          db1(i)=( z(i) * ( s(i) * z(i)  )  )

!          x(i)=( z(i) * ( s(i) * z(i) + R1 )  )

          x(i)= z(i) * s(i) * z(i) 

!+ R2 


          XMin=MIN(XMin,x(i))
          XMax=MAX(XMax,x(i))
       ENDDO

       WRITE(*,*)' NEW ',minmax(x)

44     FORMAT(A8,24(E12.6,', '))

       IF(j==1)xo_analytic=XMin
       xo_prev=xo
       xo=XMin
       !
!       sc=scale(0d0)
       sc=scale(ABS(xo))

       sc=1d0

!!$       IF( xo < ShiftSw )THEN
!!$          DO i=1,N          
!!$             x(i)=shft(x(i), dx, low_old=0d0, high_old=1d0, low_new=1d-3, high_new=1d0 )             
!!$
!!$!             x(i)=shft(x(i), dx, low_old=xo_analytic, high_old=1d0, low_new=0d0, high_new=1d0 )             
!!$          ENDDO
!!$       ENDIF

!!!       WRITE(*,44)' SHFT',X

       DO i=1,N          
         t(i) = sq_scale(x(i),sc)
!           t(i) = rugged_sq(x(i),sc,1d-3)
!           db2(i)=x(i)*t(i)**2
!          x(i) = sq_scale(x(i),sc)
       ENDDO

       WRITE(*,*)'  T  ',minmax(t)

!!$
       XMin= 1d10
       XMax=-1d10
       DO i=1,N
          XMin=MIN(XMin,x(i))
          XMax=MAX(XMax,x(i))
       ENDDO

       DO i=1,N
          R1=2d0*(RAND()-0.5d0)*SpAMM_THRESHOLD
          R2=2d0*(RAND()-0.5d0)*SpAMM_THRESHOLD
          z(i) = z(i) * t(i) + R1            
       ENDDO
       
       WRITE(*,*)'  Z  ',minmax(z)

       STOP



       IF(j>1 .AND. xo < ShiftSw )THEN
          xo_analytic=xo_analytic*(9d0/4d0)*sc
       ENDIF

       ! Convergence logic metrics
       DIFF_ONE_PREV=DIFF_ONE
       DIFF_ISQ_PREV=DIFF_ISQ
       DIFF_ONE=SUM(ABS(x(:)-1d0))/DBLE(N)
       DIFF_ISQ=SUM(ABS(z(:)-1d0/SQRT(s(:))))/DBLE(N)
       RELATIVE_DELTA=ABS(DIFF_ONE-DIFF_ONE_PREV)/DIFF_ONE
       !      
!       WRITE(*,33)j,SpAMM_THRESHOLD,XMin,scale(XMin),xo/xo_analytic,DIFF_ONE,DIFF_ISQ,RELATIVE_DELTA
31     FORMAT('it,     SpAMM,           xo,     scaling,   xo_a,         |X-I|,   |Z-1/S^.5| ')
33     FORMAT(I2,',  ',E8.2,', ',10(E12.6,', '))
       !
        IF(xo>5d-1.AND.(DIFF_ONE>DIFF_ONE_PREV.OR.RELATIVE_DELTA<5d-2))RETURN
!
    ENDDO
  END SUBROUTINE ScaledNSHalf

END module test_schulz

  program test_random_number

    USE test_schulz

    INTEGER, PARAMETER :: N=24!10**4 !876
    INTEGER :: I
    LOGICAL :: NoScale
    REAL*8, DIMENSION(1:N) :: S,X1,X2,Z1,Z2
    REAL*8 :: low1,low2,high1,high2,dx,xo,xo_prev,XMax,XMin,Lm
    !          
    CALL SetScalarSpectrum(1d10,S,N)

55 format(A5,100(E10.6,', '))
    !
    X1=S
    Z1=1

    CALL ScaledNSHalf(N,S,Z1,3d-20)
    !
    DO I=1,N
       X2(I)=Z1(I)*S(I)*Z1(I)
    ENDDO


    Z2=1
!    CALL ScaledNSHalf(N,X2,Z2,1d-4)
    !
    DELTA=0d0
    DO I=1,N
       DELTA=DELTA+ABS(1d0-(Z2(I)*Z1(I))*S(I)*(Z1(I)*Z2(I)))
    ENDDO
    WRITE(*,*)" SUM(1-Z.S.Z)/N = ",DELTA/N
    DELTA=0d0
    DO I=1,N
       DELTA=DELTA+ABS(1d0/SQRT(S(I))-Z1(I)*Z2(I))
    ENDDO
    WRITE(*,*)" SUM(1/S^.5-Z)/N = ",DELTA/N
!
  end program test_random_number

