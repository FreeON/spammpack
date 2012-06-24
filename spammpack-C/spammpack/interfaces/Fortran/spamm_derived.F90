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
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SpAMMPACK)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------

!> Defines derived types used in SpAMMPACK.
MODULE SpAMM_DERIVED

  IMPLICIT NONE

  PUBLIC

  !> Define integer of length 2.
  INTEGER, PARAMETER :: INT1 = SELECTED_INT_KIND(2)  !--Integer*1

  !> Define integer of length 2.
  INTEGER, PARAMETER :: INT2 = SELECTED_INT_KIND(4)  !--Integer*2

  !> Define integer of length 4.
  INTEGER, PARAMETER :: INT4 = SELECTED_INT_KIND(9)  !--Integer*4

  !> Define integer of length 8.
  INTEGER, PARAMETER :: INT8 = SELECTED_INT_KIND(18) !--Integer*8

  !> Define float of length 4.
  INTEGER, PARAMETER :: SpAMM_SINGLE = KIND(0.0E0)   !--Real*4

  !> Define float of length 8.
  INTEGER, PARAMETER :: SpAMM_DOUBLE = KIND(0.0D0)   !--Real*8

  !> Define the float type used for SpAMM.
#ifdef SPAMM_DOUBLE
  INTEGER, PARAMETER :: SpAMM_KIND = SpAMM_DOUBLE
#else
  INTEGER, PARAMETER :: SpAMM_KIND = SpAMM_SINGLE
#endif

  !> Define the number zero.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Zero = 0D0

  !> Define the number 1/2.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Half = 5D-1

  !> Define the number 1.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_One = 1D0

  !> Define the number 2.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Two = 2D0

  !> Define the number 4.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Four = 4D0

  !> Define the number 8.
  REAL(SpAMM_KIND), PARAMETER :: SpAMM_Eight = 8D0

  !> Quaternary tree data structure.
  TYPE SpAMM_Matrix

    !> The identifier of this object.
    INTEGER :: identifier

  END TYPE SpAMM_Matrix

END MODULE SpAMM_DERIVED
