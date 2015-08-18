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
PROGRAM SpAMM_TESTONE
  USE SpAMM_DERIVED
  USE SpAMM_GLOBALS
  USE SpAMM_MNGMENT
  USE SpAMM_ALGEBRA
  USE SpAMM_CONVERT
  USE SpAMM_PROJECT
!  USE SpAMM_PACKAGE
  IMPLICIT NONE
  INTEGER           ::  L,M,N_OLD,N,Nel,I,J,K,II,JJ,TC2_cycles,Threads
  REAL(SpAMM_KIND),ALLOCATABLE, DIMENSION(:,:) ::  A, B, C,C2,CDiff,A_NOPADDING
  REAL(SpAMM_KIND),ALLOCATABLE, DIMENSION(:,:) ::  dA, dB, dC,dP
  TYPE(QuTree),POINTER  :: qC=>NULL()
  TYPE(QuTree),POINTER  :: qP=>NULL()
  TYPE(QuTree),POINTER  :: qF=>NULL()
  TYPE(QuTree),POINTER  :: qTmp1=>NULL()
  REAL(SpAMM_KIND)          :: Occ0,TrE,RelativeErrorE,RelativeErrorN, dummy_tolerance



  REAL(DOUBLE) :: TargetTrE

  CHARACTER *100        :: Name,Buffer
  LOGICAL :: DoFilter
!---------------------------------------------------------------
  REAL(SpAMM_KIND),ALLOCATABLE, DIMENSION(:)   ::  x,p,g,h
  REAL(SpAMM_KIND)   ::  beta
!---------------------------------------------------------------



  CALL GETARG(1,NAME)

  NAME=TRIM(ADJUSTL(NAME))
  CALL GETARG(2,BUFFER)
  READ(BUFFER,*)N
  CALL GETARG(3,BUFFER)
  READ(BUFFER,*)Nel
  CALL GETARG(4,BUFFER)
  READ(BUFFER,*)TC2_cycles
  CALL GETARG(5,BUFFER)
  READ(BUFFER,*)Threads
  WRITE(*,*)' THREADS = ',THREADS
  CALL GETARG(6,BUFFER)
  READ(BUFFER,*)TargetTrE
  WRITE(*,*)NAME
  WRITE(*,*)N
  WRITE(*,*)Nel
  !--------------------------------------------------

  CALL SpAMM_Init_Globals(N,Threads)

  !
  !--------------------------------------------------
  ALLOCATE(A_NOPADDING(1:SpAMM_MATRIX_DIMENSION, &
       1:SpAMM_MATRIX_DIMENSION))
  OPEN(UNIT=66,FILE=TRIM(NAME)//".f_ortho")
  DO I=1,SpAMM_MATRIX_DIMENSION
     DO J=1,SpAMM_MATRIX_DIMENSION
        READ(66,*)L,M,A_NOPADDING(L,M)
     ENDDO
  ENDDO
  CLOSE(66)
  ALLOCATE(A(1:N,1:N))


  A=SpAMM_Zero
  A(1:SpAMM_MATRIX_DIMENSION,1:SpAMM_MATRIX_DIMENSION)= &
       A_NOPADDING(1:SpAMM_MATRIX_DIMENSION,1:SpAMM_MATRIX_DIMENSION)

  qP=>SpAMM_Convert_Dense_2_QuTree(A)
  CALL RemapSpectralBounds201(qP)
  CALL NewQuNode(qTmp1,init=.TRUE.)
  CALL AllocateFull(qTmp1)
  DO I=1,40 !TC2_cycles
     CALL SpAMM_TC2(qP,qTmp1,SpAMM_Half*FLOAT(NEl),Occ0)
     WRITE(*,*)' I ',I, Occ0*SpAMM_Two,FLOAT(NEl)
  ENDDO
  !--------------------------------------------------
  !

  CALL Multiply(qP,qF,qTmp1)
  TrE=Trace(qTmp1)
 !$OMP END SINGLE
 !$OMP END PARALLEL

  CALL SpAMM_Time_Stamp()

  RelativeErrorE=ABS(TrE-TargetTrE)/ABS(TargetTrE)
  RelativeErrorN=ABS(Occ0-SpAMM_Half*DBLE(Nel))/(SpAMM_Half*DBLE(Nel))
  !
  WRITE(*,66)N_old,N,SpAMM_BLOCK_SIZE,dummy_tolerance,RelativeErrorE,RelativeErrorN
66 FORMAT(I6,", ",I6,", ",I3,", ",E8.2,7(", ",E18.8))
  !





END PROGRAM SpAMM_TESTONE
