!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
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
!    Free Software Foundation; either version 3 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SPAMMPACK)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------

!> @brief
!! Defines conversion operation between different data structures and SpAMM.
MODULE SpAMM_CONVERT

  USE SpAMM_DERIVED
  USE SpAMM_GLOBALS
  USE SpAMM_ALGEBRA

  CONTAINS

  !=================================================================
  ! QuTree
  !=================================================================
  FUNCTION SpAMM_Convert_Dense_2_QuTree(A) RESULT(qA)

    REAL(SpAMM_KIND),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER        :: qA

    qA=>NULL()
    CALL NewQuNode(qA,init=.TRUE.)
    CALL SpAMM_Convert_Dense_2_QuTree_Recur(A,qA)
    qA%Norm=SQRT(Norm(qA))

  END FUNCTION SpAMM_Convert_Dense_2_QuTree

  RECURSIVE SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur(A,qA)
    REAL(SpAMM_KIND),DIMENSION(:,:) :: A
    TYPE(QuTree),POINTER        :: qA
    INTEGER                     :: I,J,Delta,Tier
    REAL(SpAMM_KIND)                :: Norm
    !
    I=SIZE(A,1)
    J=SIZE(A,2)
    IF(I<=SpAMM_BLOCK_SIZE.AND.J<=SpAMM_BLOCK_SIZE)THEN
       IF(I<SpAMM_BLOCK_SIZE)THEN
          WRITE(*,*)' I = ',I,' J = ',J,' Blok Siz = ',SpAMM_BLOCK_SIZE
          STOP ' LOGIC ERROR IN SpAMM: padding error '
       ELSE
!          qA%Siz=SpAMM_BLOCK_SIZE
          ALLOCATE(qA%Blok(SpAMM_BLOCK_SIZE,SpAMM_BLOCK_SIZE))
          qA%Blok(1:I,1:J)=A(1:I,1:J)
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
!!$
!!$       Tier=qA%Lev+1
!!$       qA%Quad00%Lev=Tier
!!$       qA%Quad01%Lev=Tier
!!$       qA%Quad10%Lev=Tier
!!$       qA%Quad11%Lev=Tier
!!$       !
!!$       qA%Quad00%Num=SpAMM_quadcount+1
!!$       qA%Quad01%Num=SpAMM_quadcount+2
!!$       qA%Quad10%Num=SpAMM_quadcount+3
!!$       qA%Quad11%Num=SpAMM_quadcount+4
!!$       SpAMM_quadcount=SpAMM_quadcount+4
       !
!!$       qA%Siz=(qA%Box(1,2)-qA%Box(1,1))+1
!!$       Delta=(qA%Siz-1)/2
!!$       !
!!$       qA%Quad00%Box(:,1)=qA%Box(:,1)
!!$       qA%Quad00%Box(:,2)=qA%Box(:,1)+(/Delta,Delta/)
!!$       !
!!$       qA%Quad10%Box(:,1)=(/qA%Box(1,1)+Delta+1,qA%Box(2,1)/)
!!$       qA%Quad10%Box(:,2)=(/qA%Box(1,2),qA%Box(2,1)+Delta/)
!!$       !
!!$       qA%Quad01%Box(:,1)=(/qA%Box(1,1),qA%Box(2,1)+Delta+1/)
!!$       qA%Quad01%Box(:,2)=(/qA%Box(1,1)+Delta,qA%Box(2,2)/)
!!$       !
!!$       qA%Quad11%Box(:,1)=(/qA%Box(1,1)+Delta+1,qA%Box(2,1)+Delta+1/)
!!$       qA%Quad11%Box(:,2)=qA%Box(:,2)
       !
       CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:I/2  ,1:J/2  )  , qA%Quad00 )
       CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(1:I/2  ,J/2+1:J)  , qA%Quad01 )
       CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(I/2+1:I,1:J/2  )  , qA%Quad10 )
       CALL SpAMM_Convert_Dense_2_QuTree_Recur(A(I/2+1:I,J/2+1:J)  , qA%Quad11 )
       !
    ENDIF
    !
  END SUBROUTINE SpAMM_Convert_Dense_2_QuTree_Recur

!  FUNCTION Quad2Dense(qA) RESULT(A)
!    REAL(SpAMM_KIND),DIMENSION(SpAMM_tiles*SpAMM_BLOCK_SIZE,SpAMM_tiles*SpAMM_BLOCK_SIZE) :: A
!    TYPE(QuTree),POINTER :: qA
!    A=Zero
!    CALL SpAMM_Quad2Mat(A,qA)
!  END FUNCTION Quad2Dense
!
!  RECURSIVE SUBROUTINE SpAMM_Quad2Mat(A,qA)
!    INTEGER :: I1,I2,J1,J2
!    REAL(SpAMM_KIND),DIMENSION(:,:) :: A
!    TYPE(QuTree),POINTER :: qA
!
!    IF(.NOT.ASSOCIATED(qA))RETURN
!    IF(qA%Norm==Zero)RETURN
!    !
!    IF(ALLOCATED(qA%Blok))THEN
!       A(I1:I2,J1:J2)=qA%Blok
!       RETURN
!    ELSE
!       CALL SpAMM_Quad2Mat(A,qA%Quad00)
!       CALL SpAMM_Quad2Mat(A,qA%Quad01)
!       CALL SpAMM_Quad2Mat(A,qA%Quad10)
!       CALL SpAMM_Quad2Mat(A,qA%Quad11)
!    ENDIF
!  END SUBROUTINE SpAMM_Quad2Mat

!  RECURSIVE SUBROUTINE Print_Quad(qA,PrintOnce)
!    TYPE(QuTree),POINTER :: qA
!    LOGICAL, OPTIONAL    :: PrintOnce
!
!    IF(qA%Lev>2)RETURN
!
!    IF(ASSOCIATED(qA))THEN
!       WRITE(*,111)qA%Box(:,1),qA%Box(:,2),qA%Num,qA%Siz,qA%Norm
!    ENDIF
!    IF(PRESENT(PrintOnce))THEN
!       IF(PrintOnce)RETURN
!    ENDIF
!
!    IF(ASSOCIATED(qA%Quad00))&
!         WRITE(*,111)qA%Quad00%Box(:,1),qA%Quad00%Box(:,2),qA%Quad00%Num,qA%Quad00%Siz,qA%Quad00%Norm
!    IF(ASSOCIATED(qA%Quad01))&
!         WRITE(*,111)qA%Quad01%Box(:,1),qA%Quad01%Box(:,2),qA%Quad01%Num,qA%Quad01%Siz,qA%Quad01%Norm
!    IF(ASSOCIATED(qA%Quad10))&
!         WRITE(*,111)qA%Quad10%Box(:,1),qA%Quad10%Box(:,2),qA%Quad10%Num,qA%Quad10%Siz,qA%Quad10%Norm
!    IF(ASSOCIATED(qA%Quad11))&
!         WRITE(*,111)qA%Quad11%Box(:,1),qA%Quad11%Box(:,2),qA%Quad11%Num,qA%Quad11%Siz,qA%Quad11%Norm
!    WRITE(*,*)' ========================================'
!    !
!    IF(ASSOCIATED(qA%Quad00)) &
!         CALL Print_Quad(qA%Quad00)
!    IF(ASSOCIATED(qA%Quad01)) &
!         CALL Print_Quad(qA%Quad01)
!    IF(ASSOCIATED(qA%Quad10)) &
!         CALL Print_Quad(qA%Quad10)
!    IF(ASSOCIATED(qA%Quad11)) &
!         CALL Print_Quad(qA%Quad11)
!
!111 FORMAT("Cuboid[{",I8,",",I8,":",I8,",",I8,"}, (* ",I6,", ",I6,", ",F12.6,"*)")
!
!  END SUBROUTINE Print_Quad

END MODULE SpAMM_CONVERT
