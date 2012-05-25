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
!-------------------------------------------------------------------------------
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SpAMMPack)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------

!> @mainpage
!!
!! @brief
!! The Sparse Approximate Matrix Multiply Package (SpAMMPack) is a library for
!! approximation multiplication of sparse matrices.
!!
!! @details
!! The library can be used by loading the main module SpAMM_PACKAGE which
!! includes the following sub modules:
!!
!! - SpAMM_ALGEBRA
!! - SpAMM_CONVERT
!! - SpAMM_DERIVED
!! - SpAMM_GLOBALS
!! - SpAMM_MNGMENT
!! - SpAMM_PROJECT
!!
!! Using the library from Fortran requires a
!!
!! @code
!! USE spammpack
!! @endcode
!!
!! See the documentation of @ref spammpack for more details.
!!
!! @copyright
!! This code is part of the FreeON suite of programs for linear scaling
!! electronic structure theory and ab initio molecular dynamics.
!!
!! @copyright
!! Copyright (2004). The Regents of the University of California. This
!! material was produced under U.S. Government contract W-7405-ENG-36
!! for Los Alamos National Laboratory, which is operated by the University
!! of California for the U.S. Department of Energy. The U.S. Government has
!! rights to use, reproduce, and distribute this software.  NEITHER THE
!! GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!! OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!!
!! @copyright
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by the
!! Free Software Foundation; either version 3 of the License, or (at your
!! option) any later version. Accordingly, this program is distributed in
!! the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!! PURPOSE. See the GNU General Public License at www.gnu.org for details.
!!
!! @copyright
!! While you may do as you like with this software, the GNU license requires
!! that you clearly mark derivative software.  In addition, you are encouraged
!! to return derivative works to the FreeON group for review, and possible
!! dissemination in future releases.
!!
!! @author Matt Challacombe <matt.challacombe@freeon.org>
!! @author Nicolas Bock <nicolas.bock@freeon.org>

!> The main module of SpAMMPack.
MODULE spammpack

  USE SpAMM_DERIVED
  USE SpAMM_CONVERT
  USE SpAMM_GLOBALS
  USE SpAMM_MNGMENT
  USE SpAMM_ALGEBRA
  USE SpAMM_PROJECT

CONTAINS

END MODULE spammpack
