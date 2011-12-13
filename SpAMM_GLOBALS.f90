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
MODULE SpAMM_GLOBALS
  USE SpAMM_DERIVED
  INTEGER,PARAMETER                               :: SpAMM_BLOCK_SIZE=16
#ifdef SPAMM_DOUBLE
  REAL(SpAMM_KIND),PARAMETER                      :: SpAMM_PRODUCT_TOLERANCE=1D-12
#else
  REAL(SpAMM_KIND),PARAMETER                      :: SpAMM_PRODUCT_TOLERANCE=1D-8
#endif
  REAL(SpAMM_KIND),PARAMETER                      :: SpAMM_MATRIX_TOLERANCE=1D-4*SpAMM_PRODUCT_TOLERANCE
  REAL(SpAMM_KIND),PARAMETER                      :: SpAMM_RECURSION_NORMD_CUTOFF=1E-4
  INTEGER                                         :: SpAMM_TOTAL_DEPTH
  INTEGER                                         :: SpAMM_TOTAL_LEVELS
  INTEGER                                         :: SpAMM_MATRIX_DIMENSION
  INTEGER                                         :: SpAMM_PADDED_MATRIX_DIMENSION
  INTEGER                                         :: SpAMM_THREAD_COUNT
  INTEGER                                         :: SpAMM_RECURSION_DEPTH_CUTOFF
  INTEGER,PARAMETER                               :: SpAMM_NUMBER_OF_STATS=100
  TYPE(Stats), DIMENSION(1:SpAMM_NUMBER_OF_STATS) :: SpAMM_STATS
  INTERFACE 
     FUNCTION SpAMM_IPM_GET_TIME()
       USE SpAMM_DERIVED
       REAL(SpAMM_DOUBLE) :: SpAMM_IPM_GET_TIME
     END FUNCTION SpAMM_IPM_GET_TIME
  END INTERFACE
  CONTAINS
    SUBROUTINE SpAMM_Init_Globals(N,Threads)
      IMPLICIT NONE
      INTEGER          :: N
      INTEGER,OPTIONAL :: Threads
      INTEGER          :: K,ThreadID,NThreads
      INTEGER          :: SpAMM_TILES
      !
      K=CEILING(LOG10(DBLE(N))/LOG10(2D0))
      SpAMM_MATRIX_DIMENSION=N
      SpAMM_PADDED_MATRIX_DIMENSION=2**K 
      SpAMM_TILES=CEILING(DBLE(SpAMM_PADDED_MATRIX_DIMENSION)/SpAMM_BLOCK_SIZE)
      ! Depth starts from 0:
      SpAMM_TOTAL_DEPTH=CEILING(LOG(DBLE(SpAMM_TILES))/LOG(2D0))
      ! Cutoff is tuning policy choice:
      SpAMM_RECURSION_DEPTH_CUTOFF=MAX(0,SpAMM_TOTAL_DEPTH-2)
      !
      N=SpAMM_PADDED_MATRIX_DIMENSION
!      WRITE(*,*)' SpAMM_MATRIX_DIMENSION        = ',      SpAMM_MATRIX_DIMENSION
!      WRITE(*,*)' SpAMM_PADDED_MATRIX_DIMENSION = ',SpAMM_PADDED_MATRIX_DIMENSION
#ifdef _OPENMP
      !
      IF(PRESENT(Threads))THEN
         SpAMM_THREAD_COUNT=Threads
      ELSE
         SpAMM_THREAD_COUNT=1
      ENDIF
      !
      CALL OMP_SET_DYNAMIC(.FALSE.)
      CALL OMP_SET_NESTED(.TRUE.)
      CALL OMP_SET_NUM_THREADS(SpAMM_THREAD_COUNT)
      !$OMP PARALLEL PRIVATE(ThreadID,NThreads)
      ThreadID=OMP_GET_THREAD_NUM()
      NThreads=OMP_GET_NUM_THREADS()
      WRITE(*,33)ThreadID,NThreads
      !$OMP END PARALLEL

33    FORMAT(' SpAMM_Init_Globals: Id#',I2,' checking in with ',I2,' threads ')
#endif
      SpAMM_STATS(:)%Time=0
      SpAMM_STATS(:)%Count=0
      SpAMM_STATS(:)%Routine=" "
    END SUBROUTINE SpAMM_Init_Globals

    SUBROUTINE SpAMM_Time_Stamp(Time,Routine,RoutineID)
      IMPLICIT NONE
      REAL(SpAMM_DOUBLE),OPTIONAL     :: Time
      INTEGER,OPTIONAL          :: RoutineID
      CHARACTER(LEN=*),OPTIONAL :: Routine
      REAL(SpAMM_DOUBLE)              :: SpAMM_Total_Time
      INTEGER                   :: I
      IF(.NOT.PRESENT(Time))THEN
         WRITE(*,22)
22       FORMAT(72('-'))
         DO I=1,SpAMM_NUMBER_OF_STATS
            IF(SpAMM_STATS(I)%Routine(1:1).NE." ")THEN
               WRITE(*,33)ADJUSTL(TRIM(ADJUSTL(SpAMM_STATS(I)%Routine))), &
                    SpAMM_STATS(I)%Time/DBLE(SpAMM_STATS(I)%Count),SpAMM_STATS(I)%Count
33             FORMAT(A50,' AveTime = ',F20.10,' over ',I3,' counts ') 
            ENDIF
         ENDDO
         WRITE(*,22)
         SpAMM_Total_Time=SpAMM_Zero
         DO I=1,SpAMM_NUMBER_OF_STATS
            IF(SpAMM_STATS(I)%Routine(1:1).NE." ")THEN
               SpAMM_Total_Time=SpAMM_Total_Time+SpAMM_STATS(I)%Time
               WRITE(*,34)ADJUSTL(TRIM(ADJUSTL(SpAMM_STATS(I)%Routine))), &
                    SpAMM_STATS(I)%Time
34             FORMAT(A50,' TotalTime = ',F20.10) 
            ENDIF
         ENDDO
         WRITE(*,22)
         WRITE(*,35)SpAMM_THREAD_COUNT,SpAMM_Total_Time
!         WRITE(77,35)SpAMM_THREAD_COUNT,SpAMM_Total_Time
35       FORMAT("SpAMM_SCALING ",I4,"  ", F20.10)
      ELSE
         SpAMM_STATS(RoutineID)%Time=SpAMM_STATS(RoutineID)%Time+Time
         SpAMM_STATS(RoutineID)%Count=SpAMM_STATS(RoutineID)%Count+1
         SpAMM_STATS(RoutineID)%Routine=ADJUSTL(Routine)
!         WRITE(*,44)ADJUSTL(TRIM(Routine)),Time
44       FORMAT(A50,': Time = ',F20.10,' CPUsec')
      ENDIF
    END SUBROUTINE SpAMM_Time_Stamp


  END MODULE SpAMM_GLOBALS
