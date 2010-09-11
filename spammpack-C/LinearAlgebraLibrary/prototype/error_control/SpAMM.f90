
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
  REAL(DOUBLE),PARAMETER :: Zero=0D0,Half=5D-1,One=1D0,Two=2D0
  !
  TYPE QuTree
     INTEGER                      :: Siz
     INTEGER                      :: Lev     ! Level of this box
     INTEGER                      :: Num     ! Box number
     REAL(DOUBLE)                 :: Norm
     INTEGER,      DIMENSION(2,2) :: Box
     TYPE(QuTree), POINTER        :: Quad00,Quad01,Quad10,Quad11
     REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: Blok
  END TYPE QuTree
  !
  TYPE SpAMM_cubes
     TYPE(SpAMM_cubes), POINTER  :: Next
     INTEGER                     :: Lev
     INTEGER                     :: Hash
     REAL(DOUBLE)                :: Bound
     INTEGER,     DIMENSION(3,2) :: Box
  END TYPE SpAMM_cubes
  !
  !----------------------------------------------------------------
  REAL(DOUBLE)       :: SpAMM_tolerance,SpAMM_multiplies
  INTEGER            :: SpAMM_tiles,    SpAMM_levels,     SpAMM_quadcount
  TYPE(SpAMM_cubes),POINTER :: SpAMM_stream
  INTEGER, PARAMETER :: BLOCK_SIZE=1
  !----------------------------------------------------------------
  !===============================================================================
  !  Interface blocks for generic linear algebra routines
  !===============================================================================
  INTERFACE Trace
     MODULE PROCEDURE Trace_qutree
  END INTERFACE

  INTERFACE Copy
     MODULE PROCEDURE Copy_qutree_eq_qutree
  END INTERFACE

  INTERFACE Add
     MODULE PROCEDURE Add_qutree_pls_qutree
  END INTERFACE

  INTERFACE Multiply
     MODULE PROCEDURE Multiply_qutree_tms_qutree,Multiply_qutree_tms_scalar
  END INTERFACE
  !
CONTAINS
  !
  FUNCTION Dense2Quad(A) RESULT(qA)
    REAL(DOUBLE),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER        :: qA
    !
    SpAMM_quadcount=0
    !

    CALL NewQuNode(qA,init=.TRUE.)
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
    TYPE(QuTree),POINTER :: qA
    A=Zero
    CALL SpAMM_Quad2Mat(A,qA)
  END FUNCTION Quad2Dense

  RECURSIVE SUBROUTINE SpAMM_Quad2Mat(A,qA)
    INTEGER :: I1,I2,J1,J2
    REAL(DOUBLE),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER :: qA

    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(qA%Norm==Zero)RETURN
    !
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
!    DEALLOCATE(qA)
!    NULLIFY(qA)
  END SUBROUTINE DeleteQuad

  RECURSIVE SUBROUTINE SpAMM_DeleteQuad(qA)
    TYPE(QuTree),POINTER  :: qA
    INTEGER :: Status
    IF(ALLOCATED(qA%Blok))THEN
       DEALLOCATE(qA%Blok,STAT=Status)
       NULLIFY(qA%Blok)
       SpAMM_quadcount=SpAMM_quadcount-1
    ENDIF
    IF(ASSOCIATED(qA%Quad00))THEN
       CALL SpAMM_DeleteQuad(qA%Quad00)
       DEALLOCATE(qA%Quad00)
       NULLIFY(qA%Quad00)
       SpAMM_quadcount=SpAMM_quadcount-1
    ENDIF
    IF(ASSOCIATED(qA%Quad01))THEN
       CALL SpAMM_DeleteQuad(qA%Quad01)
       DEALLOCATE(qA%Quad01)
       NULLIFY(qA%Quad01)
       SpAMM_quadcount=SpAMM_quadcount-1
    ENDIF
    IF(ASSOCIATED(qA%Quad10))THEN
       CALL SpAMM_DeleteQuad(qA%Quad10)
       DEALLOCATE(qA%Quad10)
       NULLIFY(qA%Quad10)
       SpAMM_quadcount=SpAMM_quadcount-1
    ENDIF
    IF(ASSOCIATED(qA%Quad11))THEN
       CALL SpAMM_DeleteQuad(qA%Quad11)
       DEALLOCATE(qA%Quad11)
       NULLIFY(qA%Quad11)
       SpAMM_quadcount=SpAMM_quadcount-1
    ENDIF
  END SUBROUTINE SpAMM_DeleteQuad

  SUBROUTINE NewQuNode(qA,count,init)
    LOGICAL,OPTIONAL :: init
    INTEGER,OPTIONAL :: count
    TYPE(QuTree), POINTER :: qA
    IF(PRESENT(init))THEN
       ALLOCATE(qA)
       qA%Num=1
       qA%Lev=1
       qA%Siz=SpAMM_tiles*BLOCK_SIZE
       qA%Box(:,1)=(/1,1/)
       qA%Box(:,2)=(/SpAMM_tiles,SpAMM_tiles/)*BLOCK_SIZE
    ELSE
       IF(.NOT.ASSOCIATED(qA))THEN
          ALLOCATE(qA)
       ENDIF
       IF(PRESENT(count))THEN
          count=count+1
          qA%Num=count
       ELSE
          SpAMM_quadcount=SpAMM_quadcount+1
          qA%Num=SpAMM_quadcount
       ENDIF
    ENDIF
    qA%Norm=Zero
    NULLIFY(qA%Quad00)
    NULLIFY(qA%Quad01)
    NULLIFY(qA%Quad10)
    NULLIFY(qA%Quad11)
  END SUBROUTINE NewQuNode

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

111 FORMAT("Cuboid[{",I8,",",I8,":",I8,",",I8,"}, (* ",I6,", ",I6,", ",F12.6,"*)")

  END SUBROUTINE Print_Quad

  !=================================================================
  ! QuTree trace routines: a=trace[A]
  !=================================================================
  RECURSIVE FUNCTION Trace_qutree(qA) RESULT(Trc)
    IMPLICIT NONE
    TYPE(QuTree), POINTER  :: qA
    REAL(DOUBLE) :: trc
    INTEGER :: I

!    WRITE(*,*)'TRACE , Lev =',qA%Lev,' Num = ',qA%Num,ASSOCIATED(qA%Quad00),ASSOCIATED(qA%Quad11)

    IF(qA%Siz==BLOCK_SIZE)THEN
       trc=Zero
       IF(.NOT.ALLOCATED(qA%Blok))RETURN
       DO I=1,BLOCK_SIZE
          trc=trc+qA%Blok(I,I)
       ENDDO
    ELSEIF(.NOT.(ASSOCIATED(qA%Quad00).AND.ASSOCIATED(qA%Quad11)))THEN
!       WRITE(*,*)' Z '
       Trc=Zero
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11))THEN
!       WRITE(*,*)' 00 ',ASSOCIATED(qA%Quad00)
!       WRITE(*,*)' 00 ',ASSOCIATED(qA%Quad00),ASSOCIATED(qA%Quad11)
       trc=Trace_qutree(qA%Quad00)
    ELSEIF(.NOT.ASSOCIATED(qA%Quad00))THEN
!       WRITE(*,*)' 11 '
       trc=Trace_qutree(qA%Quad11)
    ELSE
!       WRITE(*,*)' 00+11 '
       trc=Trace_qutree(qA%Quad00)+Trace_qutree(qA%Quad11)
    ENDIF
  END FUNCTION Trace_qutree
  !=================================================================
  ! QuTree filter routines: B=filter[A,tau]
  !=================================================================
  RECURSIVE SUBROUTINE Filter_qutree(qA,tau)
    IMPLICIT NONE
    TYPE(QuTree), POINTER  :: qA
    REAL(DOUBLE)           :: tau
    IF(qA%Norm<tau)THEN
       CALL DeleteQuad(qA)
    ELSE
       CALL Filter_qutree(qA%Quad00,tau)
       CALL Filter_qutree(qA%Quad01,tau)
       CALL Filter_qutree(qA%Quad10,tau)
       CALL Filter_qutree(qA%Quad11,tau)
    ENDIF
  END SUBROUTINE Filter_qutree

  !=================================================================
  ! QuTree add routines: C=A+B
  !=================================================================
  SUBROUTINE Add_qutree_pls_qutree(qA,qB,qC)
    IMPLICIT NONE
    TYPE(QuTree), POINTER       :: qA,qB,qC
    INTEGER :: add_quadcount

    IF(.NOT.ASSOCIATED(qC))THEN
       CALL NewQuNode(qC,init=.TRUE.)
    ELSE
       CALL DeleteQuad(qC)
       CALL NewQuNode(qC,init=.TRUE.)
    ENDIF
    add_quadcount=1
    CALL SpAMM_Add(qA,qB,qC,count=add_quadcount)
    qC%Norm=SQRT(qC%Norm)
  END SUBROUTINE Add_qutree_pls_qutree

  SUBROUTINE NodeMssg(Mssg,Node,count)
    CHARACTER(LEN=*) :: Mssg
    INTEGER :: count
    TYPE(QuTree), POINTER       :: Node

!    WRITE(*,11)TRIM(Mssg),Node%Lev,Node%Num,Node%Siz,count, ASSOCIATED(node), &
!     ASSOCIATED(node%Quad00),ASSOCIATED(node%Quad01),ASSOCIATED(node%Quad10),ASSOCIATED(node%Quad11)

11  FORMAT(A20,' Lev = ',I2,', Num = ',I3,', Siz = ',I4,', Cnt = ',I4,' Aloc = ',5(L2))
  END SUBROUTINE NodeMssg


  RECURSIVE SUBROUTINE SpAMM_Add(qA,qB,qC,count)
    IMPLICIT NONE
    TYPE(QuTree), POINTER    :: qA,qB,qC
    INTEGER,OPTIONAL :: count
    !

    CALL NodeMssg(' SpAMM_Add init qC ',qC,count)
    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB))THEN
       IF(.NOT.ASSOCIATED(qC))THEN
          CALL NewQuNode(qC,count=count)
          CALL NodeMssg(' SpAMM_Add newq qC ',qC,count)
       ELSE

       ENDIF
       !
       IF(qA%Siz==BLOCK_SIZE)THEN
          IF(.NOT.ALLOCATED(qC%Blok))THEN
             ALLOCATE(qC%Blok(1:BLOCK_SIZE,1:BLOCK_SIZE))
             qC%Norm=0D0
             qC%Blok=0D0
             qC%Lev=qA%Lev
             qC%Siz=BLOCK_SIZE
          ENDIF
          ! Add
          qC%Blok=qA%Blok+qB%Blok
          ! Bounds
          qC%Norm=SUM(qC%Blok(1:BLOCK_SIZE,1:BLOCK_SIZE)**2)
       ELSE
          qC%Norm=Zero

          IF(ASSOCIATED(qA%Quad00).OR.ASSOCIATED(qB%Quad00))THEN
             IF(.NOT.ASSOCIATED(qC%Quad00))THEN
                CALL NewQuNode(qC%Quad00,count=count)
             ELSE
                count=count+1
                qC%Quad00%Num=count
             ENDIF
             qC%Quad00%Lev=qA%Quad00%Lev
             qC%Quad00%Siz=qA%Quad00%Siz
             qC%Quad00%Box=qA%Quad00%Box
             CALL SpAMM_Add(qA%Quad00,qB%Quad00,qC%Quad00,count)
             qC%Norm=qC%Norm+qC%Quad00%Norm
             qC%Quad00%Norm=SQRT(qC%Quad00%Norm)
          ELSEIF(ASSOCIATED(qC%Quad00))THEN
             CALL DeleteQuad(qC%Quad00)
          ENDIF
          IF(ASSOCIATED(qA%Quad01).OR.ASSOCIATED(qB%Quad01))THEN
             IF(.NOT.ASSOCIATED(qC%Quad01))THEN
                CALL NewQuNode(qC%Quad01,count=count)
             ELSE
                count=count+1
                qC%Quad01%Num=count
             ENDIF
             qC%Quad01%Lev=qA%Quad01%Lev
             qC%Quad01%Siz=qA%Quad01%Siz
             qC%Quad01%Box=qA%Quad01%Box
             CALL SpAMM_Add(qA%Quad01,qB%Quad01,qC%Quad01,count)
             qC%Norm=qC%Norm+qC%Quad01%Norm
             qC%Quad01%Norm=SQRT(qC%Quad01%Norm)
          ELSEIF(ASSOCIATED(qC%Quad01))THEN
             CALL DeleteQuad(qC%Quad01)
          ENDIF
          IF(ASSOCIATED(qA%Quad10).OR.ASSOCIATED(qB%Quad10))THEN
             IF(.NOT.ASSOCIATED(qC%Quad10))THEN
                CALL NewQuNode(qC%Quad10,count=count)
             ELSE
                count=count+1
                qC%Quad10%Num=count
             ENDIF
             qC%Quad10%Lev=qA%Quad10%Lev
             qC%Quad10%Siz=qA%Quad10%Siz
             qC%Quad10%Box=qA%Quad10%Box
             CALL SpAMM_Add(qA%Quad10,qB%Quad10,qC%Quad10,count)
             qC%Norm=qC%Norm+qC%Quad10%Norm
             qC%Quad10%Norm=SQRT(qC%Quad10%Norm)
          ELSEIF(ASSOCIATED(qC%Quad10))THEN
             CALL DeleteQuad(qC%Quad10)
          ENDIF
          IF(ASSOCIATED(qA%Quad11).OR.ASSOCIATED(qB%Quad11))THEN
             IF(.NOT.ASSOCIATED(qC%Quad11))THEN
                CALL NewQuNode(qC%Quad11,count=count)
             ELSE
                count=count+1
                qC%Quad11%Num=count
             ENDIF
             qC%Quad11%Lev=qA%Quad11%Lev
             qC%Quad11%Siz=qA%Quad11%Siz
             qC%Quad11%Box=qA%Quad11%Box
             CALL SpAMM_Add(qA%Quad11,qB%Quad11,qC%Quad11,count)
             qC%Norm=qC%Norm+qC%Quad11%Norm
             qC%Quad11%Norm=SQRT(qC%Quad11%Norm)
          ELSEIF(ASSOCIATED(qC%Quad11))THEN
             CALL DeleteQuad(qC%Quad11)
          ENDIF
       ENDIF
    ELSEIF(.NOT.ASSOCIATED(qA))THEN
       CALL SpAMM_Copy(qB,qC,count=count)
    ELSEIF(.NOT.ASSOCIATED(qB))THEN
       CALL SpAMM_Copy(qA,qC,count=count)
    ENDIF
  END SUBROUTINE SpAMM_Add
  !=================================================================
  ! QuTree copy routines: C=A
  !=================================================================
  SUBROUTINE Copy_qutree_eq_qutree(qA,qC)
    IMPLICIT NONE
    TYPE(QuTree), POINTER       :: qA,qC
    INTEGER :: node_count
    IF(.NOT.ASSOCIATED(qC)) &
    CALL NewQuNode(qC,SpAMM_quadcount)
    qC%Num=1
    qC%Lev=1
    node_count=0
    qC%Box(:,1)=(/1,1/)
    qC%Box(:,2)=(/SpAMM_tiles,SpAMM_tiles/)*BLOCK_SIZE
    CALL SpAMM_Copy(qA,qC,count=node_count)
  END SUBROUTINE Copy_qutree_eq_qutree

  RECURSIVE SUBROUTINE SpAMM_Copy(qA,qC,count)
    IMPLICIT NONE
    TYPE(QuTree), POINTER    :: qA,qC
    INTEGER, OPTIONAL :: count
    !
    IF(.NOT.ASSOCIATED(qA))RETURN
    !
    IF(.NOT.ASSOCIATED(qC))THEN
       CALL NewQuNode(qC,count=count)
    ELSEIF(PRESENT(count))THEN
       count=count+1
       qC%Num=count
    ELSE
       SpAMM_quadcount=SpAMM_quadcount+1
       qC%Num=SpAMM_quadcount
    ENDIF
    qC%Siz=qA%Siz
    qC%Lev=qA%Lev
    qC%Box=qA%Box
    qC%Norm=qA%Norm
    CALL NodeMssg(' SpAMM_Copy qC ',qC,count)
    IF(qA%Siz==BLOCK_SIZE)THEN
       IF(.NOT.ALLOCATED(qC%Blok)) &
          ALLOCATE(qC%Blok(1:BLOCK_SIZE,1:BLOCK_SIZE))
       qC%Blok=qA%Blok
    ELSE
       IF(ASSOCIATED(qA%Quad00))THEN
          CALL SpAMM_Copy(qA%Quad00,qC%Quad00,count=count)
       ELSEIF(ASSOCIATED(qC%Quad00))THEN
          CALL DeleteQuad(qC%Quad00)
       ENDIF
       IF(ASSOCIATED(qA%Quad01))THEN
          CALL SpAMM_Copy(qA%Quad01,qC%Quad01,count=count)
       ELSEIF(ASSOCIATED(qC%Quad01))THEN
          CALL DeleteQuad(qC%Quad01)
       ENDIF
       IF(ASSOCIATED(qA%Quad10))THEN
          CALL SpAMM_Copy(qA%Quad10,qC%Quad10,count=count)
       ELSEIF(ASSOCIATED(qC%Quad10))THEN
          CALL DeleteQuad(qC%Quad10)
       ENDIF
       IF(ASSOCIATED(qA%Quad11))THEN
          CALL SpAMM_Copy(qA%Quad11,qC%Quad11,count=count)
       ELSEIF(ASSOCIATED(qC%Quad11))THEN
          CALL DeleteQuad(qC%Quad11)
       ENDIF
    ENDIF
  END SUBROUTINE SpAMM_Copy

  !=================================================================
  ! QuTree multiply routines
  !=================================================================
  ! A=A*B
  !-----------------------------------------------------------------
  RECURSIVE SUBROUTINE Multiply_qutree_tms_scalar(qA,B)
    IMPLICIT NONE
    TYPE(QuTree), POINTER    :: qA
    REAL(DOUBLE)             :: B
    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(qA%Siz==BLOCK_SIZE)THEN
       qA%Norm=qA%Norm*B
       qA%Blok=qA%Blok*B
    ELSE
       CALL Multiply_qutree_tms_scalar(qA%Quad00,B)
       CALL Multiply_qutree_tms_scalar(qA%Quad01,B)
       CALL Multiply_qutree_tms_scalar(qA%Quad10,B)
       CALL Multiply_qutree_tms_scalar(qA%Quad11,B)
       qA%Norm=qA%Norm*B
    ENDIF
  END SUBROUTINE Multiply_qutree_tms_scalar
  !-----------------------------------------------------------------
  ! C=A*B
  !-----------------------------------------------------------------
  SUBROUTINE Multiply_qutree_tms_qutree(qA,qB,qC,tolerance)
    IMPLICIT NONE
    TYPE(QuTree), POINTER       :: qA,qB,qC
    REAL(DOUBLE),OPTIONAL       :: tolerance
    REAL(DOUBLE)                :: Saved_tolerance

    CALL NodeMssg('SpAMM_mult 1 A',qA,1)
    CALL NodeMssg('SpAMM_mult 1 B',qB,1)

    IF(PRESENT(tolerance))THEN
       Saved_tolerance=SpAMM_tolerance
       SpAMM_tolerance=tolerance
    ENDIF



    IF(ASSOCIATED(qC))THEN
       CALL DeleteQuad(qC)
    ENDIF

    CALL NewQuNode(qC,init=.TRUE.)
    SpAMM_quadcount=1
    CALL SpAMM_Multiply(qC,qA,qB)
    IF(PRESENT(tolerance))THEN
       SpAMM_tolerance=Saved_tolerance
    ENDIF
  END SUBROUTINE Multiply_qutree_tms_qutree
!
  RECURSIVE SUBROUTINE SpAMM_Multiply(qC,qA,qB)
    IMPLICIT NONE
    TYPE(QuTree), POINTER :: qC,qA,qB
    REAL(DOUBLE), DIMENSION(1:BLOCK_SIZE,1:BLOCK_SIZE) :: CTmp
    LOGICAL :: do00x00,do01x10,do10x00,do11x10, &
         do00x01,do01x11,do10x01,do11x11
    ! Bounds
    qC%Siz=qA%Siz
    qC%Box(1,1)=qA%Box(2,1)
    qC%Box(1,2)=qA%Box(2,2)
    qC%Box(2,1)=qB%Box(1,1)
    qC%Box(2,2)=qB%Box(1,2)
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
       ! Bounds
       qC%Norm=qC%Norm+SUM(CTmp(1:BLOCK_SIZE,1:BLOCK_SIZE)**2)
    ELSE
       qC%Norm=Zero
       ! SpAMM=Sparse Approximate Matrix-Matrix
       do00x00=ASSOCIATED(qA%Quad00).AND.ASSOCIATED(qB%Quad00)
       do01x10=ASSOCIATED(qA%Quad01).AND.ASSOCIATED(qB%Quad10)
       do00x01=ASSOCIATED(qA%Quad00).AND.ASSOCIATED(qB%Quad01)
       do01x11=ASSOCIATED(qA%Quad01).AND.ASSOCIATED(qB%Quad11)
       do10x00=ASSOCIATED(qA%Quad10).AND.ASSOCIATED(qB%Quad00)
       do11x10=ASSOCIATED(qA%Quad11).AND.ASSOCIATED(qB%Quad10)
       do10x01=ASSOCIATED(qA%Quad10).AND.ASSOCIATED(qB%Quad01)
       do11x11=ASSOCIATED(qA%Quad11).AND.ASSOCIATED(qB%Quad11)
!       WRITE(*,11)do00x00,do01x10,do00x01,do01x11,do10x00,do11x10,do10x01,do11x11
!11     FORMAT(8(L2))
       ! SpAMM, SpAMM and SpAMM
       IF(do00x00)do00x00=do00x00.AND.(qA%Quad00%Norm*qB%Quad00%Norm>SpAMM_tolerance)
       IF(do01x10)do01x10=do01x10.AND.(qA%Quad01%Norm*qB%Quad10%Norm>SpAMM_tolerance)
       IF(do00x01)do00x01=do00x01.AND.(qA%Quad00%Norm*qB%Quad01%Norm>SpAMM_tolerance)
       IF(do01x11)do01x11=do01x11.AND.(qA%Quad01%Norm*qB%Quad11%Norm>SpAMM_tolerance)
       IF(do10x00)do10x00=do10x00.AND.(qA%Quad10%Norm*qB%Quad00%Norm>SpAMM_tolerance)
       IF(do11x10)do11x10=do11x10.AND.(qA%Quad11%Norm*qB%Quad10%Norm>SpAMM_tolerance)
       IF(do10x01)do10x01=do10x01.AND.(qA%Quad10%Norm*qB%Quad01%Norm>SpAMM_tolerance)
       IF(do11x11)do11x11=do11x11.AND.(qA%Quad11%Norm*qB%Quad11%Norm>SpAMM_tolerance)
       ! 00=00*00+01*10
       IF(do00x00.OR.do01x10)THEN
          IF(.NOT.ASSOCIATED(qC%Quad00))CALL NewQuNode(qC%Quad00,count=SpAMM_quadcount)
          IF(do00x00)CALL SpAMM_Multiply(qC%Quad00,qA%Quad00,qB%Quad00)
          IF(do01x10)CALL SpAMM_Multiply(qC%Quad00,qA%Quad01,qB%Quad10)
          qC%Norm=qC%Norm+qC%Quad00%Norm
          qC%Quad00%Norm=SQRT(qC%Quad00%Norm)
       ELSEIF(ASSOCIATED(qC%Quad00))THEN
          CALL DeleteQuad(qC%Quad00)
       ENDIF
       ! 01=00*01+01*11
       IF(do00x01.OR.do01x11)THEN
          IF(.NOT.ASSOCIATED(qC%Quad01))CALL NewQuNode(qC%Quad01,count=SpAMM_quadcount)
          IF(do00x01)CALL SpAMM_Multiply(qC%Quad01,qA%Quad00,qB%Quad01)
          IF(do01x11)CALL SpAMM_Multiply(qC%Quad01,qA%Quad01,qB%Quad11)
          qC%Norm=qC%Norm+qC%Quad01%Norm
          qC%Quad01%Norm=SQRT(qC%Quad01%Norm)
       ELSEIF(ASSOCIATED(qC%Quad01))THEN
          CALL DeleteQuad(qC%Quad01)
       ENDIF
       ! 10=10*00+11*10
       IF(do10x00.OR.do11x10)THEN
          IF(.NOT.ASSOCIATED(qC%Quad10))CALL NewQuNode(qC%Quad10,count=SpAMM_quadcount)
          IF(do10x00)CALL SpAMM_Multiply(qC%Quad10,qA%Quad10,qB%Quad00)
          IF(do11x10)CALL SpAMM_Multiply(qC%Quad10,qA%Quad11,qB%Quad10)
          qC%Norm=qC%Norm+qC%Quad10%Norm
          qC%Quad10%Norm=SQRT(qC%Quad10%Norm)
       ELSEIF(ASSOCIATED(qC%Quad10))THEN
          CALL DeleteQuad(qC%Quad10)
       ENDIF
       ! 11=10*01+11*11
       IF(do10x01.OR.do11x11)THEN
          IF(.NOT.ASSOCIATED(qC%Quad11))CALL NewQuNode(qC%Quad11,count=SpAMM_quadcount)
          IF(do10x01)CALL SpAMM_Multiply(qC%Quad11,qA%Quad10,qB%Quad01)
          IF(do11x11)CALL SpAMM_Multiply(qC%Quad11,qA%Quad11,qB%Quad11)
          qC%Norm=qC%Norm+qC%Quad11%Norm
          qC%Quad11%Norm=SQRT(qC%Quad11%Norm)
       ELSEIF(ASSOCIATED(qC%Quad11))THEN
          CALL DeleteQuad(qC%Quad11)
       ENDIF
       IF(qC%Lev==1)THEN
          qC%Norm=SQRT(qC%Norm)
       ENDIF
    ENDIF
  END SUBROUTINE SpAMM_Multiply

  SUBROUTINE TC2(P,P2,Tmp1,Norm,TrP,I)
    TYPE(QuTree),POINTER :: P,P2,Tmp1 , Tmp2
    REAL(DOUBLE) :: Norm, CR1, CR2, TrP, TrP2,    TrPd,TrPd2
    INTEGER      :: I

    REAL(DOUBLE),DIMENSION(4,4) :: dP,dP2,dTmp1, dTmp2
    !-------------------------------------------------------------------------------
    IF (I.EQ.1) THEN
       TrP=Trace(P)
    ENDIF

    CALL Multiply(P,P,P2)             ! The only multiplication is a square
    TrP2=Trace(P2)
    CR1 = ABS(TrP2-Norm)              ! CR1 = Occupation error criteria
    CR2 = ABS(2.D0*TrP - TrP2 - Norm) ! CR2 = Occupation error criteria

   WRITE(*,33)I,NORM,TrP,TrP2,CR1,CR2
33 FORMAT(I4,", ",F10.5,", ",F16.8,", ",F10.5,", ",F10.5,", ",F10.5)

   IF (CR1 < CR2) THEN               ! Too many states
      CALL Copy(P2,P)                !       P=>P2
    ELSE
       CALL Multiply(P,Two)
       CALL Multiply(P2,-One)
       CALL Add(P,P2,Tmp1)             ! P = 2P-P^2
       CALL Copy(Tmp1,P)
    ENDIF
    TrP=Trace(P)

!  CR =    2.000    2.00000000000   1.6730701243227033       0.32692987567729670       0.32692987567729759
!  CR =    2.000    1.67307012432   1.3728610345479997       0.62713896545200032       2.67207859025930805E-002
!  CR =    2.000    1.97327921409   1.8845595332763059       0.11544046672369412       6.19988949185081850E-002
!  CR =    2.000    2.06199889491   1.9985696681490941       1.43033185090590820E-003  0.12542812168792228


  END SUBROUTINE TC2

   END MODULE SpAMM_TYPES

PROGRAM SpAMM_TEST
  USE SpAMM_TYPES
  IMPLICIT NONE
  INTEGER           ::  L,M,N_OLD,N,Nel,I,J,K,II,JJ
  REAL(DOUBLE),ALLOCATABLE, DIMENSION(:,:) ::  A, B, C,C2,CDiff,A_NOPADDING
  TYPE(QuTree),POINTER  :: qC,qP,qF,qTmp1,qTmp2
  TYPE(SpAMM_cubes),POINTER  :: Stream
  REAL(DOUBLE)          :: CNorm,CDiffNorm,CErrBound,CElementErr,Opacity,Occ0,Occ1,Occ2,Occ3
  CHARACTER *100        :: Name,Buffer

  CALL GETARG(1,NAME)
  NAME=TRIM(ADJUSTL(NAME))
  CALL GETARG(2,BUFFER)
  READ(BUFFER,*)N
  CALL GETARG(3,BUFFER)
  READ(BUFFER,*)Nel
  CALL GETARG(4,BUFFER)
  READ(BUFFER,*)SpAMM_tolerance
  !
  K=CEILING(LOG10(DBLE(N))/LOG10(2D0))
  N_OLD=N
  N=2**K
  SpAMM_tiles=CEILING(DBLE(N)/BLOCK_SIZE)!*BLOCK_SIZE
  SpAMM_levels=CEILING(LOG(DBLE(SpAMM_tiles))/LOG(2D0))+1
  !
  ALLOCATE(A_NOPADDING(1:N_OLD,1:N_OLD))
  !--------------------------------------------------
  OPEN(UNIT=66,FILE=TRIM(NAME)//".f_ortho")
  DO I=1,N_OLD
     DO J=1,N_OLD
        READ(66,*)L,M,A_NOPADDING(L,M)
     ENDDO
  ENDDO
  CLOSE(66)
  ALLOCATE(A(1:N,1:N))
  A=Zero
  A(1:N_OLD,1:N_OLD)=A_NOPADDING(1:N_OLD,1:N_OLD)
  qF=>Dense2Quad(A)
  !--------------------------------------------------
  OPEN(UNIT=66,FILE=TRIM(NAME)//".p_guess_ortho")
  DO I=1,N_OLD
     DO J=1,N_OLD
        READ(66,*)L,M,A_NOPADDING(L,M)
     ENDDO
  ENDDO
  CLOSE(66)
  A=Zero
  A(1:N_OLD,1:N_OLD)=A_NOPADDING(1:N_OLD,1:N_OLD)
  qP=>Dense2Quad(A)
  !--------------------------------------------------
!  DEALLOCATE(A)
  DEALLOCATE(A_NOPADDING)
  !
  NULLIFY(qTmp1)
  NULLIFY(qTmp2)
  CALL NewQuNode(qTmp1,init=.TRUE.)
  CALL NewQuNode(qTmp2,init=.TRUE.)
  !
  Occ0 = 0.D0
  Occ1 = 0.D0
  Occ2 = 0.D0
  Occ3 = 0.D0
  DO I=1,10
    CALL TC2(qP,qTmp1,qTmp2,Half*DBLE(NEl),Occ0,I)
    Occ3 = Occ2
    Occ2 = Occ1
    Occ1 = Occ0
 ENDDO


!!$#ifdef SpAMM_print
!!$  ALLOCATE(SpAMM_stream)
!!$  Stream=>SpAMM_stream
!!$#ifdef
!!$  qC=>SpAMM(qA,qB,SpAMM_tolerance)
!!$  WRITE(44,*)'Graphics3D[{Yellow,'
!!$  DO WHILE(ASSOCIATED(Stream%Next))
!!$     Stream=>Stream%Next
!!$     IF(Stream%Lev==SpAMM_levels)THEN
!!$        Opacity=1D0
!!$     ELSE
!!$        Opacity=0.4D0*DBLE(Stream%Lev)/DBLE(SpAMM_levels)
!!$     ENDIF
!!$     WRITE(44,111)Opacity,Stream%Box(1,1),Stream%Box(2,1),Stream%Box(3,1), &
!!$          Stream%Box(1,2),Stream%Box(2,2),Stream%Box(3,2)
!!$111  FORMAT("Opacity[",F10.5,"],Cuboid[{",I6,",",I6,",",I6,"},{",I6,",",I6,",",I6,"}], ")
!!$  ENDDO
!!$
!!$  WRITE(44,*)'Boxed->False,ViewAngle -> Automatic, ViewCenter -> {0.5, 0.5, 0.5}, '
!!$  WRITE(44,*)'ViewMatrix -> Automatic, ViewPoint -> {-1.17551, -1.22153, -2.92849},'
!!$  WRITE(44,*)'ViewRange -> All, ViewVector -> Automatic,'
!!$  WRITE(44,*)'ViewVertical -> {-0.961419, -0.126619, -0.244213}]'
!!$  WRITE(*,*)' DONE SPAMM '
!!$
!!$  STOP
!!$
!!$  C2=Quad2Dense(qC)
!!$  C1=MATMUL(A,B)
!!$  CDiff=C1-C2
!!$  CNorm=SQRT(SUM(C1(1:N,1:N)**2))
!!$  CDiffNorm=SQRT(SUM(CDiff(1:N,1:N)**2))
!!$  CErrBound=SpAMM_tolerance*SQRT(SUM(A(1:N,1:N)**2))*SQRT(SUM(B(1:N,1:N)**2))
!!$
!!$  DO I=1,N
!!$     DO J=1,N
!!$        CDiff(I,J)=ABS(CDiff(I,J)) !/C1(I,J))
!!$     ENDDO
!!$  ENDDO
!!$  CElementErr=MAXVAL(CDIFF)
!!$
!!$  DO I=1,N
!!$     DO J=1,N
!!$        IF(CDiff(I,J)==CElementErr)THEN
!!$           !           WRITE(*,*)I,J
!!$           GOTO 101
!!$        ENDIF
!!$     ENDDO
!!$  ENDDO
!!$101 CONTINUE
!!$
!!$  WRITE(*,33)N,SpAMM_tolerance,CDiffNorm/CErrBound,CElementErr,SpAMM_multiplies/DBLE(SpAMM_tiles**3)
!!$33 FORMAT(I6,', ', 3(D12.6,', '),D12.6)
!!$
!!$  ! WRITE(*,*)C1(I,J),C2(I,J),(C1(I,J)-C2(I,J))/C1(I,J)
!!$
  STOP

  !
!!$  CALL DeleteQuad(qA)
!!$  CALL DeleteQuad(qB)
!!$  CALL DeleteQuad(qC)
!!$  !
!!$  DEALLOCATE(A)
!!$  DEALLOCATE(B)
!!$  DEALLOCATE(C1)
!!$  DEALLOCATE(C2)
!!$  DEALLOCATE(CDiff)

END PROGRAM SpAMM_TEST
