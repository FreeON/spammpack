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

!> Defines global parameters and variables.
MODULE SpAMM_GLOBALS
  USE SpAMM_DERIVED
  IMPLICIT NONE
  !> Define interface to the C function spamm_exit().
  INTERFACE SpAMM_Exit
    SUBROUTINE SpAMM_Exit (exitcode)
      INTEGER, INTENT(IN) :: exitcode
    END SUBROUTINE SpAMM_Exit
  END INTERFACE SpAMM_Exit
  !> Define interface to the C function spamm_trap().
  INTERFACE SpAMM_Trap
    SUBROUTINE SpAMM_Trap ()
    END SUBROUTINE SpAMM_Trap
  END INTERFACE SpAMM_Trap
  !> Define interface to the C function spamm_get_time().
  INTERFACE SpAMM_Get_Time
    FUNCTION SpAMM_Get_Time ()
      USE SpAMM_DERIVED
      REAL(SpAMM_DOUBLE) :: SpAMM_Get_Time
    END FUNCTION SpAMM_Get_Time
  END INTERFACE SpAMM_Get_Time
  !> The size of the basic submatrix blocks.
  INTEGER, PARAMETER :: SpAMM_BLOCK_SIZE = 2
  !> The SpAMM tolerance.
#ifdef SPAMM_DOUBLE
  REAL(SpAMM_KIND),PARAMETER :: SpAMM_PRODUCT_TOLERANCE = 1D-12
#else
  REAL(SpAMM_KIND),PARAMETER :: SpAMM_PRODUCT_TOLERANCE = 1D-8
#endif
  !> The "sparsification" tolerance
  REAL(SpAMM_KIND),PARAMETER :: SpAMM_MATRIX_TOLERANCE = 1D-4*SpAMM_PRODUCT_TOLERANCE
  !> The norm cutoff for tasked recursion.
  REAL(SpAMM_KIND),PARAMETER :: SpAMM_RECURSION_NORMD_CUTOFF = 1E-4
  !> The depth of the matrix tree.
  INTEGER :: SpAMM_TOTAL_DEPTH
  !> The size of the unpadded matrix.
  INTEGER :: SpAMM_MATRIX_DIMENSION
  !> The size of the padded matrix.
  INTEGER :: SpAMM_PADDED_MATRIX_DIMENSION
  !> The number of threads requested (only applicable when using OpenMP).
  INTEGER :: SpAMM_THREAD_COUNT
  !> Cutoff the tree depth at some predefined maximum depth.
  INTEGER :: SpAMM_RECURSION_DEPTH_CUTOFF
  !> Number of timer stat slots.
  INTEGER, PARAMETER :: SpAMM_NUMBER_OF_STATS = 100
  !> The timers.
  TYPE(Stats), DIMENSION(1:SpAMM_NUMBER_OF_STATS) :: SpAMM_STATS
CONTAINS

  !> Initialize global variables.
  !!
  !! @param N The matrix size. This variable gets padded and will be set to
  !! the padded matrix dimension on exit.
  !! @param Threads The number of threads. This is an optional argument and
  !! only relevant when compiled with OpenMP.
  SUBROUTINE SpAMM_Init_Globals(N, Threads)

    INTEGER, INTENT(inout) :: N
    INTEGER, OPTIONAL      :: Threads
    INTEGER                :: K
    INTEGER                :: SpAMM_TILES
#ifdef _OPENMP
    INTEGER                :: NThreads, ThreadID
#endif

    SpAMM_MATRIX_DIMENSION = N

    K=CEILING(LOG10(DBLE(N))/LOG10(2D0))

    ! Pad matrix to right size.
    SpAMM_PADDED_MATRIX_DIMENSION=2**K

    ! Number of dense matrix blocks.
    SpAMM_TILES=CEILING(DBLE(SpAMM_PADDED_MATRIX_DIMENSION)/SpAMM_BLOCK_SIZE)
    WRITE(*,*)' TILESE = ',SpAMM_TILES

    ! Depth starts from 0:
    SpAMM_TOTAL_DEPTH=CEILING(LOG(DBLE(SpAMM_TILES))/LOG(2D0))

    WRITE(*,*)' DBLE(SpAMM_PADDED_MATRIX_DIMENSION)/SpAMM_BLOCK_SIZE ',DBLE(SpAMM_PADDED_MATRIX_DIMENSION)/SpAMM_BLOCK_SIZE

    WRITE(*,*)' SPAMM_TOTAL_DEPTH = ',SPAMM_TOTAL_DEPTH
    STOP

    ! Cutoff for parallel task recursion. This is a tuning policy choice.
    SpAMM_RECURSION_DEPTH_CUTOFF=MAX(0,SpAMM_TOTAL_DEPTH-2)

    ! Set N on the way out to the padded dimension.
    N=SpAMM_PADDED_MATRIX_DIMENSION

#ifdef _OPENMP
    IF(PRESENT(Threads))THEN
      SpAMM_THREAD_COUNT=Threads
    ELSE
      SpAMM_THREAD_COUNT=1
    ENDIF

    CALL OMP_SET_DYNAMIC(.FALSE.)
    CALL OMP_SET_NESTED(.TRUE.)
    CALL OMP_SET_NUM_THREADS(SpAMM_THREAD_COUNT)

    !$OMP PARALLEL PRIVATE(ThreadID, NThreads)
    ThreadID=OMP_GET_THREAD_NUM()
    NThreads=OMP_GET_NUM_THREADS()
    !$OMP CRITICAL
    WRITE(*,33) ThreadID, NThreads
    !$OMP END CRITICAL
33  FORMAT(' SpAMM_Init_Globals: Id#',I3,' checking in with ',I3,' threads ')
    !$OMP END PARALLEL
#endif

    SpAMM_STATS(:)%Time=0
    SpAMM_STATS(:)%Count=0
    SpAMM_STATS(:)%Routine=" "

  END SUBROUTINE SpAMM_Init_Globals

  !> Reset the timer.
  SUBROUTINE SpAMM_Timer_Reset ()

    SpAMM_STATS(:)%Time=0
    SpAMM_STATS(:)%Count=0
    SpAMM_STATS(:)%Routine=" "

  END SUBROUTINE SpAMM_Timer_Reset

  !> Print out a timestamp.
  !!
  !! @detail
  !! This function keeps track of up to SpAMM_GLOBALS::SpAMM_NUMBER_OF_STATS
  !! timers. Eeach timer slot is identified by a unique RoutineID, which has to
  !! be chose by the caller in some fashion, and a string which describes the
  !! timer. Repeated calls with the same RoutineID will increment the total time
  !! in the timer slot. When called without arguments, the current state of the
  !! timers is printed.
  !!
  !! @param Time The time to print.
  !! @param Routine A string indicating the part of the code that took this
  !! amount of time.
  !! @RoutineID An ID that identifies the routine.
  SUBROUTINE SpAMM_Time_Stamp (Time, Routine, RoutineID)

    REAL(SpAMM_DOUBLE),OPTIONAL :: Time
    CHARACTER(LEN=*),OPTIONAL   :: Routine
    INTEGER,OPTIONAL            :: RoutineID

    REAL(SpAMM_DOUBLE)          :: SpAMM_Total_Time
    INTEGER                     :: I

    IF(.NOT.PRESENT(Time))THEN
      WRITE(*,22)
22    FORMAT(72('-'))
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
35    FORMAT("SpAMM_SCALING ",I4," threads ", F20.10)
    ELSE
      IF(.NOT.PRESENT(Routine) .OR. .NOT.PRESENT(RoutineID)) THEN
        WRITE(*, *) "missing Routine and/or RoutineID argument"
        CALL SpAMM_Trap()
      ENDIF

      SpAMM_STATS(RoutineID)%Time=SpAMM_STATS(RoutineID)%Time+Time
      SpAMM_STATS(RoutineID)%Count=SpAMM_STATS(RoutineID)%Count+1
      SpAMM_STATS(RoutineID)%Routine=ADJUSTL(Routine)
      !         WRITE(*,44)ADJUSTL(TRIM(Routine)),Time
      !44       FORMAT(A50,': Time = ',F20.10,' CPUsec')
    ENDIF

  END SUBROUTINE SpAMM_Time_Stamp

  !> Set the number of threads for subsequent operations.
  !!
  !! @param num_threads The number of threads to use for subsequent
  !! calculations.
  SUBROUTINE SpAMM_Set_Num_Threads (num_threads)

    INTEGER :: num_threads

    IF(num_threads >= 1) THEN
      WRITE(*, "(A,I3)") "Setting number of threads to ", num_threads
      SpAMM_THREAD_COUNT = num_threads
    ENDIF

  END SUBROUTINE SpAMM_Set_Num_Threads

END MODULE SpAMM_GLOBALS
