!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
!    LIBRARY FOR SPARSE APPROXIMATE MATRICES (libSpAMM)
!
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------
!
MODULE SpAMM_TYPES
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: INT1=SELECTED_INT_KIND(2)  !--Integer*1
  INTEGER, PARAMETER :: INT2=SELECTED_INT_KIND(4)  !--Integer*2
  INTEGER, PARAMETER :: INT4=SELECTED_INT_KIND(9)  !--Integer*4
  INTEGER, PARAMETER :: INT8=SELECTED_INT_KIND(18) !--Integer*8
  !
  INTEGER, PARAMETER :: SINGLE=KIND(0.0)           !--Real*4
  INTEGER, PARAMETER :: DOUBLE=KIND(0.D0)          !--Real*8
  REAL(DOUBLE),PARAMETER :: Zero=0D0
  !
  !----------------------------------------------------------------
  REAL(DOUBLE)       :: SpAMM_tolerance,SpAMM_multiplies
  INTEGER            :: SpAMM_tiles,    SpAMM_levels,     SpAMM_quadcount
  INTEGER, PARAMETER :: BLOCK_SIZE=2
  !----------------------------------------------------------------
  TYPE QuTree
     INTEGER                      :: Siz
     INTEGER                      :: Lev     ! Level of this box
     INTEGER                      :: Num     ! Box number
     REAL(DOUBLE)                 :: Norm
     INTEGER,      DIMENSION(2,2) :: Box
     TYPE(QuTree), POINTER        :: Quad00,Quad01,Quad10,Quad11
     REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: Blok
  END TYPE QuTree

CONTAINS

  FUNCTION Dense2Quad(A) RESULT(qA)
    REAL(DOUBLE),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER        :: qA
    !
    SpAMM_quadcount=0
    CALL NewQuNode(qA,SpAMM_quadcount)
    !
    qA%Lev=1
    qA%Num=SpAMM_quadcount
    qA%Box(:,1)=(/1,1/)
    qA%Box(:,2)=(/SpAMM_tiles,SpAMM_tiles/)*BLOCK_SIZE
    CALL SpAMM_Mat2Quad(A,qA)
    qA%Norm=SQRT(qA%Norm)
  END FUNCTION Dense2Quad

  RECURSIVE SUBROUTINE SpAMM_Mat2Quad(A,qA)
    REAL(DOUBLE),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER        :: qA
    INTEGER                     :: I,J,Delta,Tier
    REAL(DOUBLE)                :: Norm
    !
    I=SIZE(A,1)
    J=SIZE(A,2)
    IF(I<=BLOCK_SIZE.AND.J<=BLOCK_SIZE)THEN
       IF(I<BLOCK_SIZE)THEN
          STOP ' LOGIC ERROR IN SpAMM: padding error '
       ELSE
          qA%Siz=BLOCK_SIZE
          ALLOCATE(qA%Blok(BLOCK_SIZE,BLOCK_SIZE))
          qA%Blok(1:I,1:J)=A(1:I,1:J)
          qA%Norm=SUM(A(1:I,1:J)**2)
          NULLIFY(qA%Quad00)
          NULLIFY(qA%Quad01)
          NULLIFY(qA%Quad10)
          NULLIFY(qA%Quad11)
       ENDIF
       RETURN
    ELSE
       !
       ALLOCATE(qA%Quad00)
       ALLOCATE(qA%Quad01)
       ALLOCATE(qA%Quad10)
       ALLOCATE(qA%Quad11)
       !
       Tier=qA%Lev+1
       qA%Quad00%Lev=Tier
       qA%Quad01%Lev=Tier
       qA%Quad10%Lev=Tier
       qA%Quad11%Lev=Tier
       !
       qA%Quad00%Num=SpAMM_quadcount+1
       qA%Quad01%Num=SpAMM_quadcount+2
       qA%Quad10%Num=SpAMM_quadcount+3
       qA%Quad11%Num=SpAMM_quadcount+4
       SpAMM_quadcount=SpAMM_quadcount+4
       !
       qA%Siz=(qA%Box(1,2)-qA%Box(1,1))+1
       Delta=(qA%Siz-1)/2
       !
       qA%Quad00%Box(:,1)=qA%Box(:,1)
       qA%Quad00%Box(:,2)=qA%Box(:,1)+(/Delta,Delta/)
       !
       qA%Quad01%Box(:,1)=(/qA%Box(1,1)+Delta+1,qA%Box(2,1)/)
       qA%Quad01%Box(:,2)=(/qA%Box(1,2),qA%Box(2,1)+Delta/)
       !
       qA%Quad10%Box(:,1)=(/qA%Box(1,1),qA%Box(2,1)+Delta+1/)
       qA%Quad10%Box(:,2)=(/qA%Box(1,1)+Delta,qA%Box(2,2)/)
       !
       qA%Quad11%Box(:,1)=(/qA%Box(1,1)+Delta+1,qA%Box(2,1)+Delta+1/)
       qA%Quad11%Box(:,2)=qA%Box(:,2)
       !
       CALL SpAMM_Mat2Quad(A(1:I/2  ,1:J/2  )  , qA%Quad00 )
       CALL SpAMM_Mat2Quad(A(1:I/2  ,J/2+1:J)  , qA%Quad01 )
       CALL SpAMM_Mat2Quad(A(I/2+1:I,1:J/2  )  , qA%Quad10 )
       CALL SpAMM_Mat2Quad(A(I/2+1:I,J/2+1:J)  , qA%Quad11 )
       !
       qA%Norm=qA%Quad00%Norm+qA%Quad01%Norm+qA%Quad10%Norm+qA%Quad11%Norm
       !
       qA%Quad00%Norm=SQRT(qA%Quad00%Norm)
       qA%Quad10%Norm=SQRT(qA%Quad10%Norm)
       qA%Quad01%Norm=SQRT(qA%Quad01%Norm)
       qA%Quad11%Norm=SQRT(qA%Quad11%Norm)
    ENDIF
    !
  END SUBROUTINE SpAMM_Mat2Quad

  FUNCTION Quad2Dense(qA) RESULT(A)
    REAL(DOUBLE),DIMENSION(SpAMM_tiles*BLOCK_SIZE,SpAMM_tiles*BLOCK_SIZE) :: A
    TYPE(QuTree) :: qA
    CALL SpAMM_Quad2Mat(A,qA)
  END FUNCTION Quad2Dense
  RECURSIVE SUBROUTINE SpAMM_Quad2Mat(A,qA)
    INTEGER :: I1,I2,J1,J2
    REAL(DOUBLE),DIMENSION(:,:) :: A
    TYPE(QuTree) :: qA

    IF(qA%Norm==Zero)RETURN

    IF(ALLOCATED(qA%Blok))THEN
       I1=qA%Box(1,1)
       I2=qA%Box(1,2)
       J1=qA%Box(2,1)
       J2=qA%Box(2,2)
       A(I1:I2,J1:J2)=qA%Blok
       RETURN
    ELSE
       CALL SpAMM_Quad2Mat(A,qA%Quad00)
       CALL SpAMM_Quad2Mat(A,qA%Quad01)
       CALL SpAMM_Quad2Mat(A,qA%Quad10)
       CALL SpAMM_Quad2Mat(A,qA%Quad11)
    ENDIF
  END SUBROUTINE SpAMM_Quad2Mat

  SUBROUTINE DeleteQuad(qA)
    TYPE(QuTree),POINTER :: qA
    CALL SpAMM_DeleteQuad(qA)
    DEALLOCATE(qA)
  END SUBROUTINE DeleteQuad

  RECURSIVE SUBROUTINE SpAMM_DeleteQuad(qA)
    TYPE(QuTree),POINTER  :: qA
    INTEGER :: Status
    IF(ALLOCATED(qA%Blok))THEN
       DEALLOCATE(qA%Blok,STAT=Status)
       NULLIFY(qA%Blok)
       RETURN
    END IF
    IF(ASSOCIATED(qA%Quad00))THEN
       CALL SpAMM_DeleteQuad(qA%Quad00)
       DEALLOCATE(qA%Quad00)
       NULLIFY(qA%Quad00)
    ENDIF
    IF(ASSOCIATED(qA%Quad01))THEN
       CALL SpAMM_DeleteQuad(qA%Quad01)
       DEALLOCATE(qA%Quad01)
       NULLIFY(qA%Quad01)
    ENDIF
    IF(ASSOCIATED(qA%Quad10))THEN
       CALL SpAMM_DeleteQuad(qA%Quad10)
       DEALLOCATE(qA%Quad10)
       NULLIFY(qA%Quad10)
    ENDIF
    IF(ASSOCIATED(qA%Quad11))THEN
       CALL SpAMM_DeleteQuad(qA%Quad11)
       DEALLOCATE(qA%Quad11)
       NULLIFY(qA%Quad11)
    ENDIF
  END SUBROUTINE SpAMM_DeleteQuad

  RECURSIVE SUBROUTINE Print_Quad(qA)
    TYPE(QuTree),POINTER :: qA
    IF(ASSOCIATED(qA))THEN
       IF(.NOT.ASSOCIATED(qA%Quad00))RETURN
       WRITE(*,111)qA%Box(:,1),qA%Box(:,2),qA%Num,qA%Siz,qA%Norm
      ENDIF
    IF(ASSOCIATED(qA%Quad00))&
         WRITE(*,111)qA%Quad00%Box(:,1),qA%Quad00%Box(:,2),qA%Quad00%Num,qA%Quad00%Siz,qA%Quad00%Norm
    IF(ASSOCIATED(qA%Quad10))&
         WRITE(*,111)qA%Quad01%Box(:,1),qA%Quad01%Box(:,2),qA%Quad01%Num,qA%Quad01%Siz,qA%Quad10%Norm
    IF(ASSOCIATED(qA%Quad01)) &
         WRITE(*,111)qA%Quad10%Box(:,1),qA%Quad10%Box(:,2),qA%Quad10%Num,qA%Quad10%Siz,qA%Quad01%Norm
    IF(ASSOCIATED(qA%Quad11))&
         WRITE(*,111)qA%Quad11%Box(:,1),qA%Quad11%Box(:,2),qA%Quad11%Num,qA%Quad11%Siz,qA%Quad11%Norm
    WRITE(*,*)' ========================================'
    !
    IF(ASSOCIATED(qA%Quad00)) &
         CALL Print_Quad(qA%Quad00)
    IF(ASSOCIATED(qA%Quad01)) &
         CALL Print_Quad(qA%Quad01)
    IF(ASSOCIATED(qA%Quad10)) &
         CALL Print_Quad(qA%Quad10)
    IF(ASSOCIATED(qA%Quad11)) &
         CALL Print_Quad(qA%Quad11)

111 FORMAT("Cuboid[{",I6,",",I6,":",I6,",",I6,"}, (* ",I6,", ",I6,", ",F12.6,"*)")

  END SUBROUTINE Print_Quad

  FUNCTION SpAMM(qA,qB,tolerance) RESULT(qC)
    IMPLICIT NONE
    TYPE(QuTree), POINTER                :: qA,qB
    TYPE(QuTree),  POINTER      :: qC
    REAL(DOUBLE),OPTIONAL       :: tolerance
    !
    IF(PRESENT(tolerance))THEN
       SpAMM_tolerance=tolerance
    ELSE
       SpAMM_tolerance=1D-6
    ENDIF
    !
    SpAMM_quadcount=0
    CALL NewQuNode(qC,SpAMM_quadcount)
    !
    qC%Num=1
    qC%Lev=1
    SpAMM_multiplies=Zero
    qC%Box(:,1)=(/1,1/)
    qC%Box(:,2)=(/SpAMM_tiles,SpAMM_tiles/)*BLOCK_SIZE
    !
    CALL SpAMM_Multiply(qC,qA,qB)
  END FUNCTION SpAMM


  SUBROUTINE NewQuNode(qA,count)
    INTEGER :: count
    TYPE(QuTree), POINTER :: qA
    ALLOCATE(qA);  count=count+1
    qA%Num=count
    qA%Norm=Zero
    NULLIFY(qA%Quad00)
    NULLIFY(qA%Quad01)
    NULLIFY(qA%Quad10)
    NULLIFY(qA%Quad11)
  END SUBROUTINE NewQuNode

  RECURSIVE SUBROUTINE SpAMM_Multiply(qC,qA,qB)
    IMPLICIT NONE
    TYPE(QuTree), POINTER :: qC,qA,qB
    REAL(DOUBLE), DIMENSION(1:BLOCK_SIZE,1:BLOCK_SIZE) :: CTmp
    LOGICAL :: do00x00,do01x10,do10x00,do11x10, &
               do00x01,do01x11,do10x01,do11x11
    ! Bounds

    qC%Box(1,1)=qA%Box(2,1)
    qC%Box(1,2)=qA%Box(2,2)

    qC%Box(2,1)=qB%Box(1,1)
    qC%Box(2,2)=qB%Box(1,2)

!    if(qA%Box(1,1)==28.AND.QB%Box(2,2)==29)


    ! Blocks
    IF(qA%Siz==BLOCK_SIZE)THEN
       IF(.NOT.ALLOCATED(qC%Blok))THEN
          ALLOCATE(qC%Blok(1:BLOCK_SIZE,1:BLOCK_SIZE))
          qC%Norm=0D0
          qC%Blok=0D0
          qC%Lev=qA%Lev
          qC%Siz=BLOCK_SIZE
       ENDIF
       ! Count
       SpAMM_multiplies=SpAMM_multiplies+1
       ! Blocked multiply
       CTmp=MATMUL(qA%Blok,qB%Blok)
       ! Increment
       qC%Blok=qC%Blok+CTmp
!!$
!!$    IF(qA%Box(2,1)==7.AND.qB%Box(1,2)==6)THEN
!!$!    IF(qC%Box(1,1)==7.AND.qC%Box(2,2)==6)THEN
!!$!    WRITE(*,33)qA%Box,qB%Box,qC%Box,CTMP
!!$    WRITE(*,*)qC%Box(1,1),qC%Box(2,2),CTMP,qC%Blok,qA%Norm*qB%Norm
!!$ ENDIF
!!$33  format(4(I2,', '), ' || ',4(I2,', '),' || ',4(I2,', '),F20.10)
       ! Bounds
       qC%Norm=qC%Norm+SUM(CTmp(1:BLOCK_SIZE,1:BLOCK_SIZE)**2)
       RETURN
    ENDIF
    ! New nodes
    IF(.NOT.ASSOCIATED(qC%Quad00)) &
       CALL NewQuNode(qC%Quad00,SpAMM_quadcount)
    IF(.NOT.ASSOCIATED(qC%Quad01)) &
       CALL NewQuNode(qC%Quad01,SpAMM_quadcount)
    IF(.NOT.ASSOCIATED(qC%Quad10)) &
       CALL NewQuNode(qC%Quad10,SpAMM_quadcount)
    IF(.NOT.ASSOCIATED(qC%Quad11)) &
       CALL NewQuNode(qC%Quad11,SpAMM_quadcount)
    ! SpAMM=Sparse Approximate Matrix-Matrix

  !  IF(qA%Siz<=2*BLOCK_SIZE)THEN

    do00x00=qA%Quad00%Norm*qB%Quad00%Norm>SpAMM_tolerance
    do01x10=qA%Quad01%Norm*qB%Quad10%Norm>SpAMM_tolerance
    do00x01=qA%Quad00%Norm*qB%Quad01%Norm>SpAMM_tolerance
    do01x11=qA%Quad01%Norm*qB%Quad11%Norm>SpAMM_tolerance
    do10x00=qA%Quad10%Norm*qB%Quad00%Norm>SpAMM_tolerance
    do11x10=qA%Quad11%Norm*qB%Quad10%Norm>SpAMM_tolerance
    do10x01=qA%Quad10%Norm*qB%Quad01%Norm>SpAMM_tolerance
    do11x11=qA%Quad11%Norm*qB%Quad11%Norm>SpAMM_tolerance

!!$    WRITE(*,*)' ---------------------------------------------------'
!!$    WRITE(*,33)qA%Box,qB%Box,qC%Box,qA%Norm*qB%Norm
!!$    WRITE(*,*)' DO ',do00x00,do01x10,do00x01,do01x11,do10x00,do11x10,do10x01,do11x11

    ! 00=00*00+01*10
    IF(do00x00.OR.do01x10)THEN
       IF(do00x00) &
            CALL SpAMM_Multiply(qC%Quad00,qA%Quad00,qB%Quad00)
       IF(do01x10) &
            CALL SpAMM_Multiply(qC%Quad00,qA%Quad01,qB%Quad10)
    ENDIF
    ! 01=00*01+01*11
    IF(do00x01.OR.do01x11)THEN
       IF(do00x01) &
            CALL SpAMM_Multiply(qC%Quad01,qA%Quad00,qB%Quad01)
       IF(do01x11) &
            CALL SpAMM_Multiply(qC%Quad01,qA%Quad01,qB%Quad11)
    ENDIF
    ! 10=10*00+11*10
    IF(do10x00.OR.do11x10)THEN
       IF(do10x00) &
            CALL SpAMM_Multiply(qC%Quad10,qA%Quad10,qB%Quad00)
       IF(do11x10) &
            CALL SpAMM_Multiply(qC%Quad10,qA%Quad11,qB%Quad10)
    ENDIF
    ! 11=10*01+11*11
    IF(do10x01.OR.do11x11)THEN
       IF(do10x01) &
            CALL SpAMM_Multiply(qC%Quad11,qA%Quad10,qB%Quad01)
       IF(do11x11) &
            CALL SpAMM_Multiply(qC%Quad11,qA%Quad11,qB%Quad11)
    ENDIF
    ! Bounds
    qC%Norm=qC%Quad00%Norm+qC%Quad01%Norm+qC%Quad10%Norm+qC%Quad11%Norm
    qC%Quad00%Norm=SQRT(qC%Quad00%Norm)
    qC%Quad10%Norm=SQRT(qC%Quad10%Norm)
    qC%Quad01%Norm=SQRT(qC%Quad01%Norm)
    qC%Quad11%Norm=SQRT(qC%Quad11%Norm)
    IF(qC%Lev==1)THEN
       qC%Norm=SQRT(qC%Norm)
    ENDIF
    ! Done
  END SUBROUTINE SpAMM_Multiply

END MODULE SpAMM_TYPES


PROGRAM SpAMM_TEST
  USE SpAMM_TYPES
  IMPLICIT NONE
  INTEGER           ::  N,I,J,M,II,JJ
  REAL(DOUBLE),ALLOCATABLE, DIMENSION(:,:) ::  A, B, C1,C2,CDiff
  TYPE(QuTree),POINTER  :: qA,qB,qC
  REAL(DOUBLE)          :: CNorm,CDiffNorm,CErrBound,CElementErr
  CHARACTER *100        :: Buffer

  CALL GETARG(1,BUFFER)
  READ(BUFFER,*)N
  CALL GETARG(2,BUFFER)
  READ(BUFFER,*)SpAMM_tolerance
  !
  ALLOCATE(A(1:N,1:N))
  ALLOCATE(B(1:N,1:N))
  ALLOCATE(C1(1:N,1:N))
  ALLOCATE(C2(1:N,1:N))
  ALLOCATE(CDiff(1:N,1:N))

  CALL RANDOM_SEED
  CALL RANDOM_NUMBER(A)
  CALL RANDOM_NUMBER(B)

  DO I=1,N
     DO J=1,N
        !
        ! >>> COMMENTS: If you make "gamma" larger, the accounting for
        ! >>> number of multiplies becomes funny for very small thresholds
        !
        A(I,J)=EXP(-5.0D-1*ABS(I-J)) !*A(I,J)
        B(I,J)=EXP(-2.0D-1*ABS(I-J)) !*B(I,J)
     ENDDO
  ENDDO
  !
  SpAMM_tiles=CEILING(DBLE(N)/BLOCK_SIZE)!*BLOCK_SIZE
  SpAMM_levels=FLOOR(LOG(DBLE(M))/LOG(2D0))
  !
  qA=>Dense2Quad(A)
  qB=>Dense2Quad(B)
  qC=>SpAMM(qA,qB,SpAMM_tolerance)

  C2=Quad2Dense(qC)
  C1=MATMUL(A,B)
  CDiff=C1-C2
  CNorm=SQRT(SUM(C1(1:N,1:N)**2))
  CDiffNorm=SQRT(SUM(CDiff(1:N,1:N)**2))
  CErrBound=SpAMM_tolerance*SQRT(SUM(A(1:N,1:N)**2))*SQRT(SUM(B(1:N,1:N)**2))

  DO I=1,N
     DO J=1,N
        CDiff(I,J)=ABS(CDiff(I,J)) !/C1(I,J))
     ENDDO
  ENDDO
  CElementErr=MAXVAL(CDIFF)

  DO I=1,N
     DO J=1,N
        IF(CDiff(I,J)==CElementErr)THEN
           !           WRITE(*,*)I,J
           GOTO 101
        ENDIF
     ENDDO
  ENDDO
101 CONTINUE

  WRITE(*,33)N,SpAMM_tolerance,CDiffNorm/CErrBound,CElementErr,SpAMM_multiplies/DBLE(SpAMM_tiles**3)
33 FORMAT(I6,', ', 3(D12.6,', '),D12.6)

  ! WRITE(*,*)C1(I,J),C2(I,J),(C1(I,J)-C2(I,J))/C1(I,J)

  STOP

  !
  CALL DeleteQuad(qA)
  CALL DeleteQuad(qB)
  CALL DeleteQuad(qC)
  !
  DEALLOCATE(A)
  DEALLOCATE(B)
  DEALLOCATE(C1)
  DEALLOCATE(C2)
  DEALLOCATE(CDiff)

END PROGRAM SpAMM_TEST



!!$      ALLOCATE(GlblQT%Next)
!!$      ALLOCATE(GlblQT%Next%Next)
!!$      ALLOCATE(GlblQT%Next%Next%Next)
!!$      ALLOCATE(GlblQT%Next%Next%Next%Next)
!!$      !
!!$      GlblQT=>GlblQT%Next
!!$      GlblQT%Quad=>qA%Quad00
!!$      GlblQT%Hash=Hash+1*Shift
!!$      Hash00=GlblQT%Hash
!!$      !
!!$      GlblQT=>GlblQT%Next
!!$      GlblQT%Quad=>qA%Quad01
!!$      GlblQT%Hash=Hash+2*Shift
!!$      Hash01=GlblQT%Hash
!!$      !
!!$      GlblQT=>GlblQT%Next
!!$      GlblQT%Quad=>qA%Quad10
!!$      GlblQT%Hash=Hash+3*Shift
!!$      Hash10=GlblQT%Hash
!!$      !
!!$      GlblQT=>GlblQT%Next
!!$      GlblQT%Quad=>qA%Quad11
!!$      GlblQT%Hash=Hash+4*Shift
!!$      Hash11=GlblQT%Hash
!
