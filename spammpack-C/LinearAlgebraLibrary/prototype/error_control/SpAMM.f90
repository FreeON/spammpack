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
   !
   !----------------------------------------------------------------
   REAL(DOUBLE)       :: UTol
   INTEGER, PARAMETER :: BLOCK_SIZE=1
   REAL(DOUBLE)       :: NSpAMM
   !----------------------------------------------------------------

   !
   TYPE Boxx2D
      INTEGER                            :: Lev     ! Level of this box
      INTEGER                            :: Num     ! Box number
      INTEGER,        DIMENSION(2,2)     :: Box     ! Upper and lower limits of the box
   END TYPE Boxx2D

   TYPE Boxx3D
      INTEGER                            :: Lev     ! Level of this box
      INTEGER                            :: Num     ! Box number
      INTEGER,        DIMENSION(3,2)     :: Box     ! Upper and lower limits of the box
   END TYPE Boxx3D

  TYPE QuTree
     INTEGER                   :: Siz
     REAL(DOUBLE)              :: Norm
     TYPE(Boxx2D)              :: Boxx
     TYPE(QuTree),POINTER      :: Quad00,Quad01,Quad10,Quad11
     REAL(DOUBLE),DIMENSION(:,:), POINTER :: Blok
  END TYPE QuTree

  TYPE SpAMM
     REAL(DOUBLE)              :: Norm
     TYPE(Boxx3D)              :: Boxx
     TYPE(SpAMM), POINTER      :: Oct000,Oct001,Oct010,Oct100,Oct011,Oct101,Oct110,Oct111
     TYPE(QuTree),POINTER      :: A,B,C
  END TYPE SpAMM

  TYPE LQTree
     INTEGER               :: Hash
     TYPE(LQTree), POINTER :: Next
     TYPE(QuTree), POINTER :: Quad
  END TYPE LQTree

  TYPE LSpAMM
      INTEGER               :: Hash
      TYPE(LSpAMM), POINTER :: Next
      TYPE(SpAMM),  POINTER :: Octo
   END TYPE LSpAMM

  INTEGER :: Quad, Levels

  TYPE(LQTree),POINTER :: GlblQT
  TYPE(LSpAMM),POINTER :: GlblOT

  REAL(DOUBLE),PARAMETER :: Zero=0D0

  CONTAINS

    RECURSIVE SUBROUTINE SpAMM_Mat2Quad(A,qA)
      !
      REAL(DOUBLE),DIMENSION(:,:) :: A
      TYPE(QuTree),POINTER        :: qA
      INTEGER                     :: Hash,Hash00,Hash01,Hash10,Hash11,I,J,Delta,LEV,Tier,Shift
      REAL(DOUBLE)                :: Norm
      !
      I=SIZE(A,1)
      J=SIZE(A,2)
      IF(I<=BLOCK_SIZE.AND.J<=BLOCK_SIZE)THEN
         ! This is a leaf node
         ALLOCATE(qA%Blok(BLOCK_SIZE,BLOCK_SIZE))
         qA%Blok(1:I,1:J)=A(1:I,1:J)
         qA%Norm=SUM(A(1:I,1:J)**2)
         RETURN
      ENDIF
      !
      ALLOCATE(qA%Quad00)
      ALLOCATE(qA%Quad01)
      ALLOCATE(qA%Quad10)
      ALLOCATE(qA%Quad11)
      !
      Tier=qA%Boxx%Lev+1

!!$      Shift=2**(Levels-Tier)
!!$      Hash=GlblQT%Hash
      !
      qA%Quad00%Boxx%Lev=Tier
      qA%Quad01%Boxx%Lev=Tier
      qA%Quad10%Boxx%Lev=Tier
      qA%Quad11%Boxx%Lev=Tier
      !
      qA%Quad00%Boxx%Num=Quad+1
      qA%Quad01%Boxx%Num=Quad+2
      qA%Quad10%Boxx%Num=Quad+3
      qA%Quad11%Boxx%Num=Quad+4
      !
      Quad=Quad+4
      !
      qA%Siz=(qA%Boxx%Box(1,2)-qA%Boxx%Box(1,1))
      Delta=qA%Siz/2
      !
      qA%Quad00%Boxx%Box(:,1)=qA%Boxx%Box(:,1)
      qA%Quad00%Boxx%Box(:,2)=qA%Boxx%Box(:,1)+(/Delta,Delta/)
      !
      qA%Quad01%Boxx%Box(:,1)=(/qA%Boxx%Box(1,1)+Delta+1,qA%Boxx%Box(2,1)/)
      qA%Quad01%Boxx%Box(:,2)=(/qA%Boxx%Box(1,2),qA%Boxx%Box(2,1)+Delta/)
      !
      qA%Quad10%Boxx%Box(:,1)=(/qA%Boxx%Box(1,1),qA%Boxx%Box(2,1)+Delta+1/)
      qA%Quad10%Boxx%Box(:,2)=(/qA%Boxx%Box(1,1)+Delta,qA%Boxx%Box(2,2)/)
      !
      qA%Quad11%Boxx%Box(:,1)=(/qA%Boxx%Box(1,1)+Delta+1,qA%Boxx%Box(2,1)+Delta+1/)
      qA%Quad11%Boxx%Box(:,2)=qA%Boxx%Box(:,2)
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
      !
      IF(qA%Boxx%Lev==1)THEN
         qA%Norm=SQRT(qA%Norm)
      ENDIF


    END SUBROUTINE SpAMM_Mat2Quad

    RECURSIVE SUBROUTINE SpAMM_Quad2Mat(A,qA)
      INTEGER :: I1,I2,J1,J2
      REAL(DOUBLE),DIMENSION(:,:) :: A
      TYPE(QuTree) :: qA

      IF(qA%Norm==Zero)RETURN

      IF(ASSOCIATED(qA%Blok))THEN
         I1=qA%Boxx%Box(1,1)
         I2=qA%Boxx%Box(1,2)
         J1=qA%Boxx%Box(2,1)
         J2=qA%Boxx%Box(2,2)
         A(I1:I2,J1:J2)=qA%Blok
         RETURN
      ELSE
         CALL SpAMM_Quad2Mat(A,qA%Quad00)
         CALL SpAMM_Quad2Mat(A,qA%Quad01)
         CALL SpAMM_Quad2Mat(A,qA%Quad10)
         CALL SpAMM_Quad2Mat(A,qA%Quad11)
      ENDIF
    END SUBROUTINE SpAMM_Quad2Mat

    RECURSIVE SUBROUTINE Print_Quad(qA)
      TYPE(QuTree) :: qA
111   FORMAT("Cuboid[{",I6,",",I6,":",I6,",",I6,"}, (* ",I6,", ",F12.6,"*)")


      WRITE(*,111)qA%Boxx%Box(:,1),qA%Boxx%Box(:,2),qA%Boxx%Num,qA%Norm
      IF(ASSOCIATED(qA%Quad00))&
      WRITE(*,111)qA%Quad00%Boxx%Box(:,1),qA%Quad00%Boxx%Box(:,2),qA%Quad00%Boxx%Num,qA%Quad00%Norm
      IF(ASSOCIATED(qA%Quad10))&
      WRITE(*,111)qA%Quad01%Boxx%Box(:,1),qA%Quad01%Boxx%Box(:,2),qA%Quad01%Boxx%Num,qA%Quad10%Norm
      IF(ASSOCIATED(qA%Quad01)) &
      WRITE(*,111)qA%Quad10%Boxx%Box(:,1),qA%Quad10%Boxx%Box(:,2),qA%Quad10%Boxx%Num,qA%Quad01%Norm
      IF(ASSOCIATED(qA%Quad11))&
      WRITE(*,111)qA%Quad11%Boxx%Box(:,1),qA%Quad11%Boxx%Box(:,2),qA%Quad11%Boxx%Num,qA%Quad11%Norm
      WRITE(*,*)' ========================================'

      IF(qA%Siz<=BLOCK_SIZE)RETURN

      !

      CALL Print_Quad(qA%Quad00)
      IF(ASSOCIATED(qA%Quad01)) &
      CALL Print_Quad(qA%Quad01)
      IF(ASSOCIATED(qA%Quad10)) &
      CALL Print_Quad(qA%Quad10)
      IF(ASSOCIATED(qA%Quad11)) &
      CALL Print_Quad(qA%Quad11)

    END SUBROUTINE Print_Quad




END MODULE SpAMM_TYPES

PROGRAM SpAMM_TEST
  USE SpAMM_TYPES
  IMPLICIT NONE
  INTEGER,PARAMETER ::  N = 512 !2970
  INTEGER           ::  I,J,M,II,JJ
  REAL(DOUBLE),DIMENSION(N,N) ::  A, B, C1,C2,C3,A2,CDiff
  TYPE(QuTree),POINTER  :: qA,qB,qC
  TYPE(LQTree),POINTER  :: lA,lB,lC
  REAL(DOUBLE)          :: CDiffNorm
  !
  !  OPEN(UNIT=77,FILE="PP20_24056.ASCII_Xk")
  CALL RANDOM_SEED
  CALL RANDOM_NUMBER(A)
  CALL RANDOM_NUMBER(B)

  DO I=1,N
     DO J=1,N
        A(I,J)=EXP(-1.D-1*ABS(I-J))*A(I,J)
        B(I,J)=EXP(-2.0D0*ABS(I-J))*B(I,J)
        !        READ(77,*)II,JJ,A(II,JJ)
     ENDDO
  ENDDO
!!$
!!$  A(1,1)=1D0
!!$  A(1,2)=2D0
!!$  A(2,1)=3D0
!!$  A(2,2)=4D0
!!$
!!$  B(1,1)=1D0
!!$  B(1,2)=2D0
!!$  B(2,1)=3D0
!!$  B(2,2)=4D0
!!$
!!$  WRITE(*,*)A

  C1=MATMUL(A,B)
  !
  M=BLOCK_SIZE*CEILING(DBLE(N)/BLOCK_SIZE)
  Levels=FLOOR(LOG(DBLE(M))/LOG(2D0))

  WRITE(*,*)'N = ',N,'M = ',M
  !
  Quad=1
  ALLOCATE(qA)
  qA%Boxx%Num=1
  qA%Boxx%Lev=1
  qA%Boxx%Box(:,1)=(/1,1/)
  qA%Boxx%Box(:,2)=(/M,M/)
  !
  CALL SpAMM_Mat2Quad(A,qA)
!  CALL Print_Quad(qA)
!  CALL SpAMM_Quad2Mat(A2,qA)

  !
  Quad=1
  ALLOCATE(qB)
  qB%Boxx%Num=1
  qB%Boxx%Lev=1
  qB%Boxx%Box(:,1)=(/1,1/)
  qB%Boxx%Box(:,2)=(/M,M/)
  CALL SpAMM_Mat2Quad(B,qB)


  UTol=1D-6
  Quad=1
  ALLOCATE(qC)
  qC%Boxx%Num=1
  qC%Boxx%Lev=1
  qC%Boxx%Box(:,1)=(/1,1/)
  qC%Boxx%Box(:,2)=(/M,M/)
  !
  NSpAMM=Zero
  CALL SpAMM_Multiply(qC,qA,qB)
  WRITE(*,*)' Savings = ',NSpAMM/DBLE(M)**3


!  CALL Print_Quad(qC)

!!$  UTol=1D-0
!!$  Quad=1
!!$  ALLOCATE(qC)
!!$  qC%Boxx%Num=1
!!$  qC%Boxx%Lev=1
!!$  qC%Boxx%Box(:,1)=(/1,1/)
!!$  qC%Boxx%Box(:,2)=(/M,M/)
!!$  !
!!$  CALL SpAMM_Multiply(qC,qA,qB)


  !
  CALL SpAMM_Quad2Mat(C2,qC)

  CDiff=C1-C2

  CDiffNorm=SQRT(SUM(CDiff(1:BLOCK_SIZE,1:BLOCK_SIZE)**2))
  WRITE(*,*)CDiffNorm
  WRITE(*,*)UTol*SQRT(SUM(A(1:BLOCK_SIZE,1:BLOCK_SIZE)**2))*SQRT(SUM(B(1:BLOCK_SIZE,1:BLOCK_SIZE)**2))

END PROGRAM SpAMM_TEST


RECURSIVE SUBROUTINE SpAMM_Multiply(qC,qA,qB)

  USE SpAMM_TYPES
  IMPLICIT NONE
  !
  TYPE(QuTree)                :: qA,qB,qC
  !
  REAL(DOUBLE), DIMENSION(1:BLOCK_SIZE,1:BLOCK_SIZE) :: CTmp

  INTEGER :: IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2,KA1,KA2,KB1,KB2

  LOGICAL :: do00x00,do01x10,do10x00,do11x10, &
       do00x01,do01x11,do10x01,do11x11

  qC%Siz=qA%Siz
  qC%Boxx%Lev=qA%Boxx%Lev

  IA1=qA%Boxx%Box(2,1)
  IA2=qA%Boxx%Box(2,2)
  KA1=qA%Boxx%Box(1,1)
  KA2=qA%Boxx%Box(1,2)
  !A(I1:I2,J1:J2)=qA%Blok
  KB1=qB%Boxx%Box(2,1)
  KB2=qB%Boxx%Box(2,2)
  JB1=qB%Boxx%Box(1,1)
  JB2=qB%Boxx%Box(1,2)
  !B(I1:I2,J1:J2)=qA%Blok

  !    WRITE(*,*)'A =',IA1,IA2,KA1,KA2
  !    WRITE(*,*)'B =',KB1,KB2,JB1,JB2


  qC%Boxx%Box(1,1)=IA1
  qC%Boxx%Box(1,2)=IA2
  qC%Boxx%Box(2,1)=JB1
  qC%Boxx%Box(2,2)=JB2


  !          WRITE(*,*)' BoxA = ',qA%Boxx%Box
  !          WRITE(*,*)' BoxB = ',qB%Boxx%Box
  !          WRITE(*,*)' BoxC = ',qC%Boxx%Box

  IF(qA%Siz<BLOCK_SIZE)THEN
     IF(.NOT.ASSOCIATED(qC%Blok))THEN
        ALLOCATE(qC%Blok(1:BLOCK_SIZE,1:BLOCK_SIZE))
        qC%Norm=0D0
        qC%Blok=0D0
     ENDIF
     NSpAMM=NSpAMM+1D0
     CTmp=MATMUL(qA%Blok,qB%Blok)

!     IF(IA1==1.AND.JB1==1)THEN

     qC%Blok=qC%Blok+CTmp
     qC%Norm=qC%Norm+SUM(CTmp(1:BLOCK_SIZE,1:BLOCK_SIZE)**2)

!     WRITE(*,*)' k = ',KA1,IA1,IA2,JB1,JB2,qA%Norm,qB%Norm,qC%Blok,qC%Norm


!  ENDIF
     RETURN
  ENDIF
  !
  IF(.NOT.ASSOCIATED(qC%Quad00))THEN
     ALLOCATE(qC%Quad00)
     qC%Quad00%Norm=Zero
  ENDIF
  IF(.NOT.ASSOCIATED(qC%Quad01))THEN
     ALLOCATE(qC%Quad01)
     qC%Quad01%Norm=Zero
  ENDIF
  IF(.NOT.ASSOCIATED(qC%Quad10))THEN
     ALLOCATE(qC%Quad10)
     qC%Quad10%Norm=Zero
  ENDIF
  IF(.NOT.ASSOCIATED(qC%Quad11))THEN
     ALLOCATE(qC%Quad11)
     qC%Quad11%Norm=Zero
  ENDIF
!
!!$  do00x00=.TRUE.
!!$  do01x10=.TRUE.
!!$  do00x01=.TRUE.
!!$  do01x11=.TRUE.
!!$  do10x00=.TRUE.
!!$  do11x10=.TRUE.
!!$  do10x01=.TRUE.
!!$  do11x11=.TRUE.


  do00x00=qA%Quad00%Norm*qB%Quad00%Norm>UTol
  do01x10=qA%Quad01%Norm*qB%Quad10%Norm>UTol
  do00x01=qA%Quad00%Norm*qB%Quad01%Norm>UTol
  do01x11=qA%Quad01%Norm*qB%Quad11%Norm>UTol
  do10x00=qA%Quad10%Norm*qB%Quad00%Norm>UTol
  do11x10=qA%Quad11%Norm*qB%Quad10%Norm>UTol
  do10x01=qA%Quad10%Norm*qB%Quad01%Norm>UTol
  do11x11=qA%Quad11%Norm*qB%Quad11%Norm>UTol

  ! 00=00*00+01*10
  IF(do00x00.OR.do01x10)THEN
     IF(do00x00)&
          CALL SpAMM_Multiply(qC%Quad00,qA%Quad00,qB%Quad00)
     IF(do01x10)&
          CALL SpAMM_Multiply(qC%Quad00,qA%Quad01,qB%Quad10)
  ENDIF
  ! 01=00*01+01*11
  IF(do00x01.OR.do01x11)THEN
     IF(do00x01) &
          CALL SpAMM_Multiply(qC%Quad01,qA%Quad01,qB%Quad11)
     IF(do01x11) &
          CALL SpAMM_Multiply(qC%Quad01,qA%Quad00,qB%Quad01)
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
  qC%Norm=qC%Quad00%Norm+qC%Quad01%Norm+qC%Quad10%Norm+qC%Quad11%Norm

  qC%Quad00%Norm=SQRT(qC%Quad00%Norm)
  qC%Quad10%Norm=SQRT(qC%Quad10%Norm)
  qC%Quad01%Norm=SQRT(qC%Quad01%Norm)
  qC%Quad11%Norm=SQRT(qC%Quad11%Norm)
  !
  IF(qC%Boxx%Lev==1)THEN
     qC%Norm=SQRT(qC%Norm)
  ENDIF

END SUBROUTINE SpAMM_Multiply



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
