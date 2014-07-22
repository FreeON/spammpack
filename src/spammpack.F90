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
!! @author Nicolas Bock nicolas.bock@freeon.org

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
