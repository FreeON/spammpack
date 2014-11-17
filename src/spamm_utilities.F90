!> Defines utilitiy functions.
!!
!! @copyright
!!
!! Copyright (c) 2014, Los Alamos National Laboratory
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice, this
!! list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!! this list of conditions and the following disclaimer in the documentation
!! and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its contributors
!! may be used to endorse or promote products derived from this software without
!! specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! @author Matt Challacombe matt.challacombe@freeon.org
!! @author Nicolas Bock nicolasbock@freeon.org
module spamm_utilities

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  !> Fatal error.
  integer, parameter :: FATAL = -1

  !> Debug messages.
  integer, parameter :: DEBUG = 2

  !> @deprecated Number of timer stat slots.
  INTEGER, PARAMETER :: SpAMM_NUMBER_OF_STATS = 100

  !> A type for performance measurements.
  TYPE Stats

    !> The time.
    real(kind(0d0)) :: Time

    !> Some count.
    INTEGER :: Count

    !> The name of a function.
    CHARACTER(LEN=50) :: Routine

  END TYPE Stats

  !> @deprecated The timers.
  TYPE(Stats), DIMENSION(SpAMM_NUMBER_OF_STATS) :: SpAMM_STATS

  !> Interface to to_string functions.
  interface to_string
    module procedure int_to_string
    module procedure double_array_to_string
    module procedure single_to_string
    module procedure double_to_string
    module procedure bitree_to_string
    module procedure qutree_to_string
    module procedure spamm_matrix_order_1_to_string
    module procedure spamm_matrix_2nd_order_to_string
  end interface to_string

contains

  !> Convert an integer to a string.
  !!
  !! @param i The integer
  !!
  !! @return The string representation.
  function int_to_string (i) result(str_i)

    integer, intent(in) :: i

    character(len = 100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_i
#else
    character(len = 100) :: str_i
#endif

    write(temp, *) i
    str_i = trim(adjustl(temp))

  end function int_to_string

  !> Convert a single precision real to a string.
  !!
  !! @param x The real
  !!
  !! @return The string representation.
  function single_to_string (x) result(str_x)

    real(kind(0e0)), intent(in) :: x
    character(len = 100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_x
#else
    character(len = 100), allocatable :: str_x
#endif

    write(temp, "(ES16.8E3)") x
    str_x = trim(adjustl(temp))

  end function single_to_string

  !> Convert an array of double precision reals to a string.
  !!
  !! @param x The real array.
  !!
  !! @return The string representation.
  function double_array_to_string (x) result(str_x)

    real(kind(0d0)), intent(in) :: x(:)
    character(len = 1000) :: temp
    character(len = 100) :: format_string
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_x
#else
    character(len = 100), allocatable :: str_x
#endif
    integer :: i

    write(format_string, "(A,I4,A)") "(", size(x), "ES16.8E3)"
    write(temp, format_string) (x(i), i = 1, size(x))
    str_x = trim(adjustl(temp))

  end function double_array_to_string

  !> Convert a double precision real to a string.
  !!
  !! @param x The real
  !!
  !! @return The string representation.
  function double_to_string (x) result(str_x)

    real(kind(0d0)), intent(in) :: x
    character(len = 100) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_x
#else
    character(len = 100), allocatable :: str_x
#endif

    write(temp, "(ES16.8E3)") x
    str_x = trim(adjustl(temp))

  end function double_to_string

  !> Convert a bitree node to a string.
  !!
  !! @param q The bitree node.
  !!
  !! @return The string representation.
  function bitree_to_string (q) result(str_q)

    use spamm_types

    type(bitree), pointer, intent(in) :: q
    character(len = 200) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_q
#else
    character(len = 200), allocatable :: str_q
#endif

    if(associated(q)) then
      write(temp, "(A)") "b["// &
        to_string(q%i_lower)//":"// &
        to_string(q%i_upper)//"]"

      write(temp, "(A)") trim(temp)//", norm = "//to_string(q%norm)

      if(allocated(q%vect)) then
        write(temp, "(A)") trim(temp)//", vect["// &
          int_to_string(size(q%vect, 1))//"]"
      else
        write(temp, "(A)") trim(temp)//", vect not allocated"
      endif

      if(associated(q%sect1)) then
        write(temp, "(A)") trim(temp)//", sect1"
      endif

      if(associated(q%sect2)) then
        write(temp, "(A)") trim(temp)//", sect2"
      endif

    else
      write(temp, "(A)") "q not associated"
    endif

    str_q = trim(adjustl(temp))

  end function bitree_to_string

  !> Convert a qutree node to a string.
  !!
  !! @param q The qutree node.
  !!
  !! @return The string representation.
  function qutree_to_string (q) result(str_q)

    use spamm_types

    type(qutree), pointer, intent(in) :: q
    character(len = 200) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_q
#else
    character(len = 200), allocatable :: str_q
#endif

    if(associated(q)) then
      write(temp, "(A)") "q["// &
        to_string(q%i_lower)//":"// &
        to_string(q%i_upper)//","// &
        to_string(q%j_lower)//":"// &
        to_string(q%j_upper)//"]"

      write(temp, "(A)") trim(temp)//", norm = "//to_string(q%norm)

      if(allocated(q%blok)) then
        write(temp, "(A)") trim(temp)//", blok["// &
          int_to_string(size(q%blok, 1))//","// &
          int_to_string(size(q%blok, 2))//"]"
      else
        write(temp, "(A)") trim(temp)//", block not allocated"
      endif

      if(associated(q%quad11)) then
        write(temp, "(A)") trim(temp)//", quad11"
      endif

      if(associated(q%quad12)) then
        write(temp, "(A)") trim(temp)//", quad12"
      endif

      if(associated(q%quad21)) then
        write(temp, "(A)") trim(temp)//", quad21"
      endif

      if(associated(q%quad22)) then
        write(temp, "(A)") trim(temp)//", quad22"
      endif

    else
      write(temp, "(A)") "q not associated"
    endif

    str_q = trim(adjustl(temp))

  end function qutree_to_string

  !> Convert a spamm_matrix_order_1 to a string.
  !!
  !! @param V The vector.
  !!
  !! @return The string representation.
  function spamm_matrix_order_1_to_string (V) result(str_q)

    use spamm_types

    type(spamm_matrix_order_1), pointer, intent(in) :: V
    character(len = 200) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_q
#else
    character(len = 200), allocatable :: str_q
#endif

    if(associated(V)) then
      write(temp, "(A)") "N = "//to_string(V%N)
      write(temp, "(A)") trim(temp)//", N_padded = "//to_string(V%N_padded)
      write(temp, "(A)") trim(temp)//", depth = "//to_string(V%depth)
      write(temp, "(A)") trim(temp)//", norm = "//to_string(V%norm)
    else
      write(temp, "(A)") "V not associated"
    endif

    str_q = trim(adjustl(temp))

  end function spamm_matrix_order_1_to_string

  !> Convert a spamm_matrix_2nd_order to a string.
  !!
  !! @param A The matrix.
  !!
  !! @return The string representation.
  function spamm_matrix_2nd_order_to_string (A) result(str_q)

    use spamm_types

    type(spamm_matrix_2nd_order), pointer, intent(in) :: A
    character(len = 200) :: temp
#ifdef HAVE_DEFERRED_STRING_LENGTH
    character(len = :), allocatable :: str_q
#else
    character(len = 200), allocatable :: str_q
#endif

    if(associated(A)) then
      write(temp, "(A)") "M = "//to_string(A%M)
      write(temp, "(A)") trim(temp)//", N = "//to_string(A%N)
      write(temp, "(A)") trim(temp)//", N_padded = "//to_string(A%N_padded)
      write(temp, "(A)") trim(temp)//", depth = "//to_string(A%depth)
      write(temp, "(A)") trim(temp)//", norm = "//to_string(A%norm)
    else
      write(temp, "(A)") "A not associated"
    endif

    str_q = trim(adjustl(temp))

  end function spamm_matrix_2nd_order_to_string

  !> Print a log message.
  !!
  !! @param level The message level. The global debug_level has to be greater or equal to the level.
  !! @param message The message to print.
  subroutine write_log (level, message)

    integer, intent(in) :: level
    character(len = *), intent(in) :: message

    if(SPAMM_DEBUG_LEVEL >= level) then
      if(level < 0) then
        write(*, "(A)") "[FATAL] "//message
        error stop
      else
        write(*, "(A)") message
      endif
    endif

  end subroutine write_log

  !> The timer.
  !! @return The time passed since some point in time.
  real(kind(0d0)) function spamm_get_time ()

  end function spamm_get_time

  !> @deprecated Initialize global variables.
  !!
  !! @param N The matrix size. This variable gets padded and will be set to
  !! the padded matrix dimension on exit.
  !! @param Threads The number of threads. This is an optional argument and
  !! only relevant when compiled with OpenMP.
  !> Reset the timer.
  SUBROUTINE SpAMM_Timer_Reset ()

    SpAMM_STATS(:)%Time=0
    SpAMM_STATS(:)%Count=0
    SpAMM_STATS(:)%Routine=" "

  END SUBROUTINE SpAMM_Timer_Reset

  !> Print out a timestamp.
  !!
  !! This function keeps track of up to spamm_types::spamm_number_of_stats
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

    REAL(kind(0d0)), OPTIONAL :: Time
    CHARACTER(LEN=*), OPTIONAL :: Routine
    INTEGER, OPTIONAL :: RoutineID

    REAL(kind(0d0)) :: SpAMM_Total_Time
    INTEGER :: I

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
      SpAMM_Total_Time = 0
      DO I=1,SpAMM_NUMBER_OF_STATS
        IF(SpAMM_STATS(I)%Routine(1:1).NE." ")THEN
          SpAMM_Total_Time=SpAMM_Total_Time+SpAMM_STATS(I)%Time
          WRITE(*,34)ADJUSTL(TRIM(ADJUSTL(SpAMM_STATS(I)%Routine))), &
            SpAMM_STATS(I)%Time
          34             FORMAT(A50,' TotalTime = ',F20.10)
        ENDIF
      ENDDO
      WRITE(*,22)
#ifdef _OPENMP
      WRITE(*,35) omp_get_num_threads(), SpAMM_Total_Time
#else
      WRITE(*,35) 1, SpAMM_Total_Time
#endif
35    FORMAT("SpAMM ",I4," threads ", F20.10)
    ELSE
      IF(.NOT.PRESENT(Routine) .OR. .NOT.PRESENT(RoutineID)) THEN
        WRITE(*, *) "missing Routine and/or RoutineID argument"
        error stop
      ENDIF

      SpAMM_STATS(RoutineID)%Time=SpAMM_STATS(RoutineID)%Time+Time
      SpAMM_STATS(RoutineID)%Count=SpAMM_STATS(RoutineID)%Count+1
      SpAMM_STATS(RoutineID)%Routine=ADJUSTL(Routine)
      !         WRITE(*,44)ADJUSTL(TRIM(Routine)),Time
      !44       FORMAT(A50,': Time = ',F20.10,' CPUsec')
    ENDIF

  END SUBROUTINE SpAMM_Time_Stamp

end module spamm_utilities
