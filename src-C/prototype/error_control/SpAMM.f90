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
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SPAMMPACK)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------
MODULE SpAMM_TYPES
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: INT1=SELECTED_INT_KIND(2)  !--Integer*1
  INTEGER, PARAMETER :: INT2=SELECTED_INT_KIND(4)  !--Integer*2
  INTEGER, PARAMETER :: INT4=SELECTED_INT_KIND(9)  !--Integer*4
  INTEGER, PARAMETER :: INT8=SELECTED_INT_KIND(18) !--Integer*8
  !
  INTEGER, PARAMETER :: SINGLE=KIND(0.0)           !--Real*4
  INTEGER, PARAMETER :: DOUBLE=KIND(0.0D0)         !--Real*8
  REAL(DOUBLE),PARAMETER :: Zero=0D0,Half=5D-1,One=1D0,Two=2D0,Eight=8D0
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
  REAL(DOUBLE)       :: SpAMM_tolerance, SpAMM_multiplies, SpAMM_average_multiplies, &
                        SpAMM_dimension
  INTEGER            :: SpAMM_tiles, SpAMM_levels, SpAMM_quadcount
  TYPE(SpAMM_cubes),POINTER :: SpAMM_stream
!  INTEGER, PARAMETER :: SpAMM_BLOCK_SIZE=4
  INTEGER :: SpAMM_BLOCK_SIZE

  REAL(DOUBLE) :: GlobalDropped,GlobalError

  !----------------------------------------------------------------
  !===============================================================================
  !  Interface blocks for generic linear algebra routines
  !===============================================================================
  INTERFACE Trace
     MODULE PROCEDURE Trace_qutree, Trace_dns
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
    qA%Norm=SQRT(SpAMM_NormReduce(qA))
    CALL Filter_qutree(qA,1D-16)
  END FUNCTION Dense2Quad

  RECURSIVE SUBROUTINE SpAMM_Mat2Quad(A,qA)
    REAL(DOUBLE),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER        :: qA
    INTEGER                     :: I,J,Delta,Tier
    REAL(DOUBLE)                :: Norm
    !
    I=SIZE(A,1)
    J=SIZE(A,2)
    IF(I<=SpAMM_BLOCK_SIZE.AND.J<=SpAMM_BLOCK_SIZE)THEN
       IF(I<SpAMM_BLOCK_SIZE)THEN
          WRITE(*,*)' I = ',I,' J = ',J,' Blok Siz = ',SpAMM_BLOCK_SIZE
          STOP ' LOGIC ERROR IN SpAMM: padding error '
       ELSE
          qA%Siz=SpAMM_BLOCK_SIZE
          ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))
          qA%Blok(1:I,1:J)=A(1:I,1:J)
!          qA%Norm=SUM(A(1:I,1:J)**2)
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
       qA%Quad10%Box(:,1)=(/qA%Box(1,1)+Delta+1,qA%Box(2,1)/)
       qA%Quad10%Box(:,2)=(/qA%Box(1,2),qA%Box(2,1)+Delta/)
       !
       qA%Quad01%Box(:,1)=(/qA%Box(1,1),qA%Box(2,1)+Delta+1/)
       qA%Quad01%Box(:,2)=(/qA%Box(1,1)+Delta,qA%Box(2,2)/)
       !
       qA%Quad11%Box(:,1)=(/qA%Box(1,1)+Delta+1,qA%Box(2,1)+Delta+1/)
       qA%Quad11%Box(:,2)=qA%Box(:,2)
       !
       CALL SpAMM_Mat2Quad(A(1:I/2  ,1:J/2  )  , qA%Quad00 )
       CALL SpAMM_Mat2Quad(A(1:I/2  ,J/2+1:J)  , qA%Quad01 )
       CALL SpAMM_Mat2Quad(A(I/2+1:I,1:J/2  )  , qA%Quad10 )
       CALL SpAMM_Mat2Quad(A(I/2+1:I,J/2+1:J)  , qA%Quad11 )
!       !
!       qA%Norm=qA%Quad00%Norm+qA%Quad01%Norm+qA%Quad10%Norm+qA%Quad11%Norm
!       !
!       qA%Quad00%Norm=SQRT(qA%Quad00%Norm)
!       qA%Quad01%Norm=SQRT(qA%Quad01%Norm)
!       qA%Quad10%Norm=SQRT(qA%Quad10%Norm)
!       qA%Quad11%Norm=SQRT(qA%Quad11%Norm)
    ENDIF
    !
  END SUBROUTINE SpAMM_Mat2Quad

  FUNCTION Quad2Dense(qA) RESULT(A)
    REAL(DOUBLE),DIMENSION(SpAMM_tiles*SpAMM_BLOCK_SIZE,SpAMM_tiles*SpAMM_BLOCK_SIZE) :: A
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
    DEALLOCATE(qA)
    NULLIFY(qA)
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
       qA%Siz=SpAMM_tiles*SpAMM_BLOCK_SIZE
       qA%Box(:,1)=(/1,1/)
       qA%Box(:,2)=(/SpAMM_tiles,SpAMM_tiles/)*SpAMM_BLOCK_SIZE
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
       WRITE(*,111)qA%Box(:,1),qA%Box(:,2),qA%Num,qA%Siz,qA%Norm
    ENDIF

    IF(ASSOCIATED(qA%Quad00))&
         WRITE(*,111)qA%Quad00%Box(:,1),qA%Quad00%Box(:,2),qA%Quad00%Num,qA%Quad00%Siz,qA%Quad00%Norm
    IF(ASSOCIATED(qA%Quad01))&
         WRITE(*,111)qA%Quad01%Box(:,1),qA%Quad01%Box(:,2),qA%Quad01%Num,qA%Quad01%Siz,qA%Quad01%Norm
    IF(ASSOCIATED(qA%Quad10))&
         WRITE(*,111)qA%Quad10%Box(:,1),qA%Quad10%Box(:,2),qA%Quad10%Num,qA%Quad10%Siz,qA%Quad10%Norm
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

    IF(qA%Siz==SpAMM_BLOCK_SIZE)THEN
       trc=Zero
       IF(.NOT.ALLOCATED(qA%Blok))RETURN
       DO I=1,SpAMM_BLOCK_SIZE
          trc=trc+qA%Blok(I,I)
       ENDDO
    ELSEIF(.NOT.ASSOCIATED(qA%Quad00).AND..NOT.ASSOCIATED(qA%Quad11))THEN
       Trc=Zero
   ELSEIF(.NOT.ASSOCIATED(qA%Quad11))THEN
       trc=Trace_qutree(qA%Quad00)
    ELSEIF(.NOT.ASSOCIATED(qA%Quad00))THEN
       trc=Trace_qutree(qA%Quad11)
    ELSE
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
    IF(.NOT.ASSOCIATED(qA))RETURN
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
    REAL(DOUBLE) :: norm
    INTEGER :: add_quadcount
    IF(.NOT.ASSOCIATED(qC))THEN
       CALL NewQuNode(qC,init=.TRUE.)
    ELSE
       CALL DeleteQuad(qC)
       CALL NewQuNode(qC,init=.TRUE.)
    ENDIF
    add_quadcount=1
    CALL SpAMM_Add(qA,qB,qC,count=add_quadcount)
    qC%Norm=SQRT(SpAMM_NormReduce(qC))
    CALL Filter_qutree(qC,1D-16)
  END SUBROUTINE Add_qutree_pls_qutree

  RECURSIVE SUBROUTINE SpAMM_Add(qA,qB,qC,count)
    IMPLICIT NONE
    TYPE(QuTree), POINTER    :: qA,qB,qC
    INTEGER,OPTIONAL :: count
    !
    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB))THEN
       IF(.NOT.ASSOCIATED(qC))THEN
          CALL NewQuNode(qC,count=count)
       ELSE
          count=count+1
          qC%Num=count
       ENDIF
       qC%Lev=qA%Lev
       qC%Siz=qA%Siz
       qC%Box=qA%Box
       !
       IF(qA%Siz==SpAMM_BLOCK_SIZE)THEN
          IF(.NOT.ALLOCATED(qC%Blok))THEN
             ALLOCATE(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
             qC%Blok=0D0
          ENDIF
          ! Add
          IF(ALLOCATED(qA%Blok).AND.ALLOCATED(qB%Blok))THEN
             qC%Blok=qA%Blok+qB%Blok
             qC%Lev=qA%Lev
          ELSEIF(ALLOCATED(qA%Blok))THEN
             qC%Blok=qA%Blok
             qC%Lev=qA%Lev
          ELSEIF(ALLOCATED(qB%Blok))THEN
             qC%Blok=qB%Blok
             qC%Lev=qB%Lev
          ELSE
             CALL DeleteQuad(qC)
          ENDIF
       ELSE
          CALL SpAMM_Add(qA%Quad00,qB%Quad00,qC%Quad00,count)
          CALL SpAMM_Add(qA%Quad01,qB%Quad01,qC%Quad01,count)
          CALL SpAMM_Add(qA%Quad10,qB%Quad10,qC%Quad10,count)
          CALL SpAMM_Add(qA%Quad11,qB%Quad11,qC%Quad11,count)
       ENDIF
    ELSEIF(.NOT.ASSOCIATED(qA))THEN
       CALL SpAMM_Copy(qB,qC,count=count)
    ELSEIF(.NOT.ASSOCIATED(qB))THEN
       CALL SpAMM_Copy(qA,qC,count=count)
    ELSE
       CALL DeleteQuad(qC)
    ENDIF
  END SUBROUTINE SpAMM_Add
  !=================================================================
  ! QuTree copy routines: C=A
  !=================================================================
  FUNCTION Copy_qutree_eq_qutree(qA) RESULT(qC)
    IMPLICIT NONE
    TYPE(QuTree), POINTER       :: qA,qC
    INTEGER :: node_count
    CALL NewQuNode(qC,init=.TRUE.)
    qC%Num=1
    qC%Lev=1
    node_count=1
    qC%Box(:,1)=(/1,1/)
    qC%Box(:,2)=(/SpAMM_tiles,SpAMM_tiles/)*SpAMM_BLOCK_SIZE
    CALL SpAMM_Copy(qA,qC,count=node_count)
  END FUNCTION Copy_qutree_eq_qutree

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
    IF(qA%Siz==SpAMM_BLOCK_SIZE.AND.ALLOCATED(qA%Blok))THEN
       IF(.NOT.ALLOCATED(qC%Blok)) &
          ALLOCATE(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
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
    IF(qA%Siz==SpAMM_BLOCK_SIZE.AND.ALLOCATED(qA%Blok))THEN
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
    REAL(DOUBLE)                :: Saved_tolerance,norm,AbsError,RelError,TargetError,ABNorm
    INTEGER                     :: SpAMM_count
    REAL(DOUBLE),ALLOCATABLE,DIMENSION(:,:) :: dA,dB,dC

    IF(PRESENT(tolerance))THEN
       Saved_tolerance=SpAMM_tolerance
       SpAMM_tolerance=tolerance
    ENDIF
    CALL DeleteQuad(qC)
    CALL NewQuNode(qC,init=.TRUE.)
    SpAMM_count=1
    SpAMM_multiplies=0

    GlobalDropped=Zero
    GlobalError=Zero


    CALL SpAMM_Multiply2(qC,qA,qB,count=SpAMM_count)
    qC%Norm=SQRT(SpAMM_NormReduce(qC))
    CALL Filter_qutree(qc,1D-16)


!!$    ABNorm=qA%Norm*qB%Norm
!!$    TargetError=SpAMM_tolerance/SQRT(Eight)
!!$    CALL SpAMM_Multiply(qC,qA,qB,TargetError,count=SpAMM_count)
!!$    qC%Norm=SQRT(SpAMM_NormReduce(qC))
!!$    CALL Filter_qutree(qc,1D-16)
!!$    ALLOCATE(dA(1:SpAMM_tiles*SpAMM_BLOCK_SIZE,1:SpAMM_tiles*SpAMM_BLOCK_SIZE))
!!$    ALLOCATE(dB(1:SpAMM_tiles*SpAMM_BLOCK_SIZE,1:SpAMM_tiles*SpAMM_BLOCK_SIZE))
!!$    ALLOCATE(dC(1:SpAMM_tiles*SpAMM_BLOCK_SIZE,1:SpAMM_tiles*SpAMM_BLOCK_SIZE))
!!$
!!$    dA=Quad2Dense(qA)
!!$    dB=Quad2Dense(qB)
!!$    dC=MATMUL(dA,dB)
!!$    AbsError=NormReduce_dns(dC-Quad2Dense(qC))
!!$    RelError=AbsError/ABNorm !NormReduce_dns(dC)
!!$    DEALLOCATE(dA)
!!$    DEALLOCATE(dB)
!!$    DEALLOCATE(dC)
!!$!    WRITE(*,*)' True  Error          = ',TrueError
!!$    WRITE(*,33)SpAMM_tolerance,GlobalDropped,SpAMM_multiplies,AbsError,RelError
!!$33  FORMAT('Tol = ',D12.6,', Drops = ',D12.6,', Mults = ',D12.6,', Abs Error = ',D12.6,', Rel Error = ',D12.6)

    !
    IF(PRESENT(tolerance))THEN
       SpAMM_tolerance=Saved_tolerance
    ENDIF
  END SUBROUTINE Multiply_qutree_tms_qutree
  !
  RECURSIVE FUNCTION SpAMM_NormReduce(qC) RESULT(norm)
    IMPLICIT NONE
    TYPE(QuTree), POINTER :: qC
    REAL(DOUBLE) :: norm
    ! Associated
    IF(.NOT.ASSOCIATED(qC))THEN
       norm=Zero
       RETURN
    ELSEIF(qC%Siz==SpAMM_BLOCK_SIZE.AND.ALLOCATED(qC%Blok))THEN
       norm=SUM(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE)**2)
       qC%Norm=SQRT(norm)
    ELSE
       norm=SpAMM_NormReduce(qC%Quad00) + &
            SpAMM_NormReduce(qC%Quad01) + &
            SpAMM_NormReduce(qC%Quad10) + &
            SpAMM_NormReduce(qC%Quad11)
       qC%Norm=SQRT(norm)
    ENDIF
  END FUNCTION SpAMM_NormReduce

  FUNCTION NormReduce_dns(A) RESULT(norm)
    REAL(DOUBLE), DIMENSION(:,:) :: A
    REAL(DOUBLE) :: norm
    INTEGER :: I,J
    norm=zero
    DO I=1,SIZE(A,1)
       DO J=1,SIZE(A,2)
          norm=norm+A(I,J)**2
       ENDDO
    ENDDO
    norm=SQRT(norm)
  END FUNCTION NormReduce_dns

  RECURSIVE SUBROUTINE SpAMM_Multiply(qC,qA,qB,TotalError,count)
    IMPLICIT NONE
    TYPE(QuTree), POINTER :: qC,qA,qB
    INTEGER     :: count,LocalCount,LocalDropped
    REAL(DOUBLE) :: TotalError,LocalError,ErrorEstimate
    REAL(DOUBLE):: norm00x00,norm00x01,norm01x11, &
                   norm10x00,norm10x01,norm11x11,norm01x10,norm11x10
    REAL(DOUBLE):: frac00x00,frac01x10,frac00x01,frac01x11, &
                   frac10x00,frac11x10,frac10x01,frac11x11
    LOGICAL     :: do00x00,do01x10,do10x00,do11x10, &
                   do00x01,do01x11,do10x01,do11x11
    ! Bounds
    qC%Lev=qA%Lev
    qC%Siz=qA%Siz
    qC%Box(1,:)=qA%Box(1,:)
    qC%Box(2,:)=qB%Box(2,:)
    ! Blocks
    IF(qA%Siz==SpAMM_BLOCK_SIZE)THEN
       ! Allocate
       IF(.NOT.ALLOCATED(qC%Blok))THEN
          ALLOCATE(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
          qC%Blok=Zero
       END IF
       ! Count
       SpAMM_multiplies=SpAMM_multiplies+1
       ! Accumulate
       qC%Blok=qC%Blok+MATMUL(qA%Blok,qB%Blok)
    ELSE
       ! Nodes exist?
       do00x00=ASSOCIATED(qA%Quad00).AND.ASSOCIATED(qB%Quad00)
       do01x10=ASSOCIATED(qA%Quad01).AND.ASSOCIATED(qB%Quad10)
       do00x01=ASSOCIATED(qA%Quad00).AND.ASSOCIATED(qB%Quad01)
       do01x11=ASSOCIATED(qA%Quad01).AND.ASSOCIATED(qB%Quad11)
       do10x00=ASSOCIATED(qA%Quad10).AND.ASSOCIATED(qB%Quad00)
       do11x10=ASSOCIATED(qA%Quad11).AND.ASSOCIATED(qB%Quad10)
       do10x01=ASSOCIATED(qA%Quad10).AND.ASSOCIATED(qB%Quad01)
       do11x11=ASSOCIATED(qA%Quad11).AND.ASSOCIATED(qB%Quad11)
       !
       norm00x00=Zero
       norm01x10=Zero
       norm00x01=Zero
       norm01x11=Zero
       norm10x00=Zero
       norm11x10=Zero
       norm10x01=Zero
       norm11x11=Zero
       !
       IF(do00x00)norm00x00=qA%Quad00%Norm*qB%Quad00%Norm
       IF(do01x10)norm01x10=qA%Quad01%Norm*qB%Quad10%Norm
       IF(do00x01)norm00x01=qA%Quad00%Norm*qB%Quad01%Norm
       IF(do01x11)norm01x11=qA%Quad01%Norm*qB%Quad11%Norm
       IF(do10x00)norm10x00=qA%Quad10%Norm*qB%Quad00%Norm
       IF(do11x10)norm11x10=qA%Quad11%Norm*qB%Quad10%Norm
       IF(do10x01)norm10x01=qA%Quad10%Norm*qB%Quad01%Norm
       IF(do11x11)norm11x11=qA%Quad11%Norm*qB%Quad11%Norm
       ! SpAMM criteria
       IF(do00x00)do00x00=do00x00.AND.norm00x00>TotalError
       IF(do01x10)do01x10=do01x10.AND.norm01x10>TotalError
       IF(do00x01)do00x01=do00x01.AND.norm00x01>TotalError
       IF(do01x11)do01x11=do01x11.AND.norm01x11>TotalError
       IF(do10x00)do10x00=do10x00.AND.norm10x00>TotalError
       IF(do11x10)do11x10=do11x10.AND.norm11x10>TotalError
       IF(do10x01)do10x01=do10x01.AND.norm10x01>TotalError
       IF(do11x11)do11x11=do11x11.AND.norm11x11>TotalError
       ! Counting the next generation
       LocalCount=0
       IF(do00x00)LocalCount=LocalCount+1
       IF(do01x10)LocalCount=LocalCount+1
       IF(do00x01)LocalCount=LocalCount+1
       IF(do01x11)LocalCount=LocalCount+1
       IF(do10x00)LocalCount=LocalCount+1
       IF(do11x10)LocalCount=LocalCount+1
       IF(do10x01)LocalCount=LocalCount+1
       IF(do11x11)LocalCount=LocalCount+1
       ! Check for done by division
       IF(LocalCount==0)RETURN
       ! global debug/development counters
       GlobalDropped=GlobalDropped+(8-LocalCount)
       GlobalError=GlobalError+LocalError
       ! Statistical estimate of error partitioning
       TotalError=TotalError/SQRT(DBLE(LocalCount))
       ! 00=00*00+01*10
       IF(do00x00.OR.do01x10)THEN
          IF(.NOT.ASSOCIATED(qC%Quad00))CALL NewQuNode(qC%Quad00,count=count)
          IF(do00x00)THEN
             frac00x00=TotalError
             CALL SpAMM_Multiply(qC%Quad00,qA%Quad00,qB%Quad00,frac00x00,count)
          ENDIF
          IF(do01x10)THEN
             frac01x10=TotalError
             CALL SpAMM_Multiply(qC%Quad00,qA%Quad01,qB%Quad10,frac01x10,count)
          ENDIF
       ENDIF
       ! 01=00*01+01*11
       IF(do00x01.OR.do01x11)THEN
          IF(.NOT.ASSOCIATED(qC%Quad01))CALL NewQuNode(qC%Quad01,count=count)
          IF(do00x01)THEN
             frac00x01=TotalError
             CALL SpAMM_Multiply(qC%Quad01,qA%Quad00,qB%Quad01,frac00x01,count)
          ENDIF
          IF(do01x11)THEN
             frac01x11=TotalError
             CALL SpAMM_Multiply(qC%Quad01,qA%Quad01,qB%Quad11,frac01x11,count)
          ENDIF
       ENDIF
       ! 10=10*00+11*10
       IF(do10x00.OR.do11x10)THEN
          IF(.NOT.ASSOCIATED(qC%Quad10))CALL NewQuNode(qC%Quad10,count=count)
          IF(do10x00)THEN
             frac10x00=TotalError
             CALL SpAMM_Multiply(qC%Quad10,qA%Quad10,qB%Quad00,frac10x00,count)
          ENDIF
          IF(do11x10)THEN
             frac11x10=TotalError
             CALL SpAMM_Multiply(qC%Quad10,qA%Quad11,qB%Quad10,frac11x10,count)
          ENDIF
       ENDIF
       ! 11=10*01+11*11
       IF(do10x01.OR.do11x11)THEN
          IF(.NOT.ASSOCIATED(qC%Quad11))CALL NewQuNode(qC%Quad11,count=count)
          IF(do10x01)THEN
             frac10x01=TotalError
             CALL SpAMM_Multiply(qC%Quad11,qA%Quad10,qB%Quad01,frac10x01,count)
          ENDIF
          IF(do11x11)THEN
             frac11x11=TotalError
             CALL SpAMM_Multiply(qC%Quad11,qA%Quad11,qB%Quad11,frac11x11,count)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE SpAMM_Multiply

  RECURSIVE SUBROUTINE SpAMM_Multiply2(qC,qA,qB,count)
    IMPLICIT NONE
    TYPE(QuTree), POINTER :: qC,qA,qB
    INTEGER :: count
    REAL(DOUBLE) :: norm
    REAL(DOUBLE), DIMENSION(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE) :: CTmp
    ! Associated
    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB))THEN
    ! Bounds
    IF(qA%Norm*qB%Norm<SpAMM_tolerance)RETURN
    ! Allocate
    IF(.NOT.ASSOCIATED(qC))CALL NewQuNode(qC,count=count)
    ! Boxing
     qC%Siz=qA%Siz
     qC%Lev=qA%Lev
     qC%Box(1,:)=qA%Box(1,:)
     qC%Box(2,:)=qB%Box(2,:)
       ! Blocks
       IF(qA%Siz==SpAMM_BLOCK_SIZE)THEN
          ! Allocate
          IF(.NOT.ALLOCATED(qC%Blok))THEN
             ALLOCATE(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
             qC%Blok=Zero
          END IF
          ! Count
          SpAMM_multiplies=SpAMM_multiplies+1
          ! Accumulate
          qC%Blok=qC%Blok+MATMUL(qA%Blok,qB%Blok)
       ELSE
          ! 00=00*00+01*10
          CALL SpAMM_Multiply2(qC%Quad00,qA%Quad00,qB%Quad00,count=count)
          CALL SpAMM_Multiply2(qC%Quad00,qA%Quad01,qB%Quad10,count=count)
          ! 01=00*01+01*11
          CALL SpAMM_Multiply2(qC%Quad01,qA%Quad00,qB%Quad01,count=count)
          CALL SpAMM_Multiply2(qC%Quad01,qA%Quad01,qB%Quad11,count=count)
          ! 10=10*00+11*10
          CALL SpAMM_Multiply2(qC%Quad10,qA%Quad10,qB%Quad00,count=count)
          CALL SpAMM_Multiply2(qC%Quad10,qA%Quad11,qB%Quad10,count=count)
          ! 11=10*01+11*11
          CALL SpAMM_Multiply2(qC%Quad11,qA%Quad10,qB%Quad01,count=count)
          CALL SpAMM_Multiply2(qC%Quad11,qA%Quad11,qB%Quad11,count=count)
       ENDIF
    ENDIF
  END SUBROUTINE SpAMM_Multiply2

  FUNCTION Trace_dns(A) RESULT(trc)
    REAL(DOUBLE) :: trc
    REAL(DOUBLE),DIMENSION(:,:) :: A
    INTEGER I
    trc=Zero
    DO I=1,SIZE(A,1)
       trc=trc+A(I,I)
    ENDDO
  END FUNCTION Trace_dns

  SUBROUTINE TC2_test(P,P2,Tmp1,Norm,TrP,I)
    TYPE(QuTree),POINTER :: P,P2,Tmp1 , Tmp2
    REAL(DOUBLE) :: Norm, CR1, CR2, TrP, TrP2,    TrdP,TrdP2, sum
    INTEGER      :: I,  II,J,KK

    INTEGER, PARAMETER:: NN=1024
    REAL(DOUBLE),DIMENSION(NN,NN) :: dP,dP2,dTmp1, dTmp2,dP2b
    !-------------------------------------------------------------------------------

    dP=Quad2Dense(P)

    IF (I.EQ.1) THEN
       TrP=Trace(P)
       TrdP=Trace(dP)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !       TrP=TrdP
       IF(ABS(TrP-TrdP)>SpAMM_tolerance*1D1)THEN
          WRITE(*,*)' 1TrP = ',TrP,TrdP
          STOP
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ENDIF
    !
    CALL Multiply(P,P,P2)             ! The only multiplication is a square
    dP2=MATMUL(dP,dP)


!    P2=>Dense2Quad(Quad2Dense(P2))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dP2b=Quad2Dense(P2)
    DO II=1,NN
       DO J=1,NN
          IF(ABS(dP2(II,J)-dP2b(II,J))>SpAMM_tolerance*1D1)THEN
             WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
             WRITE(*,*)' sparse P*P'

             sum=Zero
             DO KK=1,NN
                IF(ABS(dP(II,KK)*dP(KK,J))>1D-1*SpAMM_tolerance)THEN
                   WRITE(*,*)dP(II,KK)*dP(KK,J),sum
                   sum=sum+dP(II,KK)*dP(KK,J)
                ENDIF
             ENDDO
             WRITE(*,*)' SpAMM_tolerance = ',SpAMM_tolerance
            WRITE(*,*)'ARAW',II,J,dP2(II,J),dP2b(II,J),dP2(II,J)-dP2b(II,J)
            WRITE(*,*)'ASUM',II,J,sum,dP2b(II,J),dP2(II,J)-sum

!             CALL Print_quad(P2)
 !            P=>Dense2Quad(dP2)
             WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
             WRITE(*,*)' dense P*P'


 !            CALL Print_quad(P)
             STOP
          ENDIF
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    TrP2=Trace(P2)
    TrdP2=Trace(dP2)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !    TrP2=TrdP2
    IF(ABS(TrP2-TrdP2)>SpAMM_tolerance)THEN
       WRITE(*,*)'TrP = ',TrP2,TrdP2
       STOP "P2 STOP "
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    CR1 = ABS(TrP2-Norm)              ! CR1 = Occupation error criteria
    CR2 = ABS(2.D0*TrP - TrP2 - Norm) ! CR2 = Occupation error criteria

    WRITE(*,33)I,NORM,TrP,SpAMM_multiplies/SpAMM_dimension**3
33  FORMAT(I4,", N =",F12.5,", Tr(P)=",F12.5,", O(N)/N^3 = ",F8.5)

    IF (CR1 < CR2) THEN               ! Too many states
      P=>Copy(P2)
      dP=dP2
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dP2b=Quad2Dense(P)
       DO II=1,NN
          DO J=1,NN
             IF(ABS(dP(II,J)-dP2b(II,J))>SpAMM_tolerance*1D1)THEN
                WRITE(*,*)' B ',II,J,dP(II,J)-dP2b(II,J)
                STOP
             ENDIF
          ENDDO
       ENDDO
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ELSE
       CALL Multiply(P,Two)
       CALL Multiply(P2,-One)
       CALL Add(P,P2,Tmp1)             ! P = 2P-P^2
       P=>Copy(Tmp1)
       dP=Two*dP-dP2
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       dP2b=Quad2Dense(Tmp1)
       DO II=1,NN
          DO J=1,NN
             IF(ABS(dP(II,J)-dP2b(II,J))>SpAMM_tolerance*1D1)THEN
                WRITE(*,*)II,J,dP(II,J)-dP2b(II,J)
                STOP
             ENDIF
          ENDDO
       ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ENDIF
    TrP=Trace(P)
    TrdP=Trace(dP)

!    TrP=TrdP

    IF(ABS(TrP-TrdP)>SpAMM_tolerance*1D1)THEN
       WRITE(*,*)'TrP = ',TrP,TrdP
       STOP
    ENDIF

  END SUBROUTINE TC2_TEST

  SUBROUTINE TC2(P,P2,Tmp1,Norm,TrP,I)
    TYPE(QuTree),POINTER :: P,P2,Tmp1 , Tmp2
    REAL(DOUBLE) :: Norm, CR1, CR2, TrP, TrP2,    TrdP,TrdP2
    INTEGER      :: I,  II,J
    !-------------------------------------------------------------------------------
    IF (I.EQ.1)TrP=Trace(P)
    CALL Multiply(P,P,P2)             ! The only multiplication is a square
    TrP2=Trace(P2)
    CR1 = ABS(TrP2-Norm)              ! CR1 = Occupation error criteria
    CR2 = ABS(2.D0*TrP - TrP2 - Norm) ! CR2 = Occupation error criteria
!    WRITE(*,33)I,NORM,TrP,SpAMM_multiplies/SpAMM_tiles**3
33  FORMAT(I4,", N =",F18.8,", Tr(P)=",F18.8,", O(N)/N^3 = ",F8.5)
    IF (CR1 < CR2) THEN               ! Too many states
      P=>Copy(P2)
    ELSE
       CALL Multiply(P,Two)
       CALL Multiply(P2,-One)
       CALL Add(P,P2,Tmp1)             ! P = 2P-P^2
       P=>Copy(Tmp1)
    ENDIF
    TrP=Trace(P)
!    CALL Filter_qutree(P,SpAMM_tolerance*1D-2)
    !
  END SUBROUTINE TC2


  SUBROUTINE TC2_filter(P,P2,Tmp1,Norm,TrP,I)
    TYPE(QuTree),POINTER :: P,P2,Tmp1 , Tmp2
    REAL(DOUBLE) :: Norm, CR1, CR2, TrP, TrP2,    TrdP,TrdP2
    INTEGER      :: I,  II,J
    !-------------------------------------------------------------------------------
    IF(I.EQ.1)THEN
       CALL Filter_qutree(P,SpAMM_tolerance)
       TrP=Trace(P)
    ENDIF
    CALL Multiply_qutree_tms_qutree(P,P,P2,tolerance=1D-100)
    TrP2=Trace(P2)
    CR1 = ABS(TrP2-Norm)              ! CR1 = Occupation error criteria
    CR2 = ABS(2.D0*TrP - TrP2 - Norm) ! CR2 = Occupation error criteria
!    WRITE(*,33)I,NORM,TrP,SpAMM_multiplies/SpAMM_tiles**3
33  FORMAT(I4,", N =",F12.8,", Tr(P)=",F12.8,", O(N)/N^3 = ",F8.5)
    IF (CR1 < CR2) THEN               ! Too many states
       P=>Copy(P2)
    ELSE
       CALL Multiply(P,Two)
       CALL Multiply(P2,-One)
       CALL Add(P,P2,Tmp1)             ! P = 2P-P^2
       P=>Copy(Tmp1)
    ENDIF
    CALL Filter_qutree(P,SpAMM_tolerance)
    TrP=Trace(P)
  END SUBROUTINE TC2_FILTER

END MODULE SpAMM_TYPES

PROGRAM SpAMM_TEST
  USE SpAMM_TYPES
  IMPLICIT NONE
  INTEGER           ::  L,M,N_OLD,N,Nel,I,J,K,II,JJ,TC2_cycles
  REAL(DOUBLE),ALLOCATABLE, DIMENSION(:,:) ::  A, B, C,C2,CDiff,A_NOPADDING
  TYPE(QuTree),POINTER  :: qC,qP,qF,qTmp1,qTmp2
  TYPE(SpAMM_cubes),POINTER  :: Stream
  REAL(DOUBLE)          :: CNorm,CDiffNorm,CErrBound,CElementErr,Opacity,Occ0,Occ1,Occ2,Occ3,TrE, &
                           RelativeErrorE,RelativeErrorN,TargetTrE
  CHARACTER *100        :: Name,Buffer
  LOGICAL :: DoFilter

  CALL GETARG(1,NAME)
  NAME=TRIM(ADJUSTL(NAME))
  CALL GETARG(2,BUFFER)
  READ(BUFFER,*)N
  CALL GETARG(3,BUFFER)
  READ(BUFFER,*)Nel
  CALL GETARG(4,BUFFER)
  READ(BUFFER,*)TC2_cycles
  CALL GETARG(5,BUFFER)
  READ(BUFFER,*)SpAMM_tolerance
  CALL GETARG(6,BUFFER)
  READ(BUFFER,*)DoFilter
  CALL GETARG(7,BUFFER)
  READ(BUFFER,*)TargetTrE
  !--------------------------------------------------
  K=CEILING(LOG10(DBLE(N))/LOG10(2D0))
  N_OLD=N
  N=2**K
  SpAMM_dimension=N

  SpAMM_BLOCK_SIZE=4
!!!  SpAMM_BLOCK_SIZE=N



  SpAMM_tiles=CEILING(DBLE(N)/SpAMM_BLOCK_SIZE)
  SpAMM_levels=CEILING(LOG(DBLE(SpAMM_tiles))/LOG(2D0))+1
  ! ~ ~
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
  !  CALL Print_quad(qP)
  !--------------------------------------------------
  DEALLOCATE(A)
  DEALLOCATE(A_NOPADDING)
  !
  NULLIFY(qTmp1)
  NULLIFY(qTmp2)
  CALL NewQuNode(qTmp1,init=.TRUE.)
  CALL NewQuNode(qTmp2,init=.TRUE.)
  !--------------------------------------------------
  SpAMM_average_multiplies=Zero
  IF(.NOT.DoFilter)THEN
!     SpAMM_tolerance=SpAMM_tolerance*DBLE(Nel)/Two
  ENDIF
  !--------------------------------------------------
  Occ0 = 0.D0;Occ1 = 0.D0;Occ2 = 0.D0;Occ3 = 0.D0
  DO I=1,TC2_cycles
     IF(DoFilter)THEN
        CALL TC2_filter(qP,qTmp1,qTmp2,Half*FLOAT(NEl),Occ0,I)
     ELSE
        CALL TC2(qP,qTmp1,qTmp2,Half*FLOAT(NEl),Occ0,I)
     ENDIF
     Occ3 = Occ2;Occ2 = Occ1;Occ1 = Occ0
     SpAMM_average_multiplies=SpAMM_average_multiplies+SpAMM_multiplies
  ENDDO

!  SpAMM_average_multiplies=SpAMM_average_multiplies/DBLE(TC2_cycles)
!  WRITE(*,*)' AV MULTIPLIES = ',SpAMM_average_multiplies

  !


  !
  IF(DoFilter)THEN
     CALL Filter_qutree(qF,SpAMM_tolerance)
     CALL Multiply_qutree_tms_qutree(qP,qF,qTmp1,tolerance=1D-100)
  ELSE
!     CALL Multiply_qutree_tms_qutree(qP,qF,qTmp1,tolerance=1D-100)
     CALL Multiply(qP,qF,qTmp1)
  ENDIF
  TrE=Trace(qTmp1)
  RelativeErrorE=ABS(TrE-TargetTrE)/ABS(TargetTrE)
  RelativeErrorN=ABS(Occ0-Half*DBLE(Nel))/(Half*DBLE(Nel))

  SpAMM_average_multiplies=SpAMM_average_multiplies+SpAMM_multiplies
  SpAMM_average_multiplies=SpAMM_average_multiplies/DBLE(TC2_cycles+1)


!!!  WRITE(*,*)"TrE = ",TrE

  !
  WRITE(*,66)N_old,N,SpAMM_BLOCK_SIZE,SpAMM_tolerance,RelativeErrorE,RelativeErrorN, &
             SpAMM_average_multiplies
66 FORMAT(I6,", ",I6,", ",I3,", ",E8.2,7(", ",E18.8))
  !
END PROGRAM SpAMM_TEST
