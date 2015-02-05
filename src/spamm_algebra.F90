!> Defines operations on SpAMM trees.
!!
!! @copyright
!!
!! Copyright (c) 2015, Los Alamos National Laboratory
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
module spamm_algebra

  use spamm_globals
  use spamm_management
  use spamm_real_precision
  use spamm_types
  use spamm_utilities

#ifdef _OPENMP
  use omp_lib
#endif

!  INCLUDE "spamm_utility_macros.h"

  implicit none

  PRIVATE

  PUBLIC :: Multiply
  PUBLIC :: Trace
  PUBLIC :: Add
  PUBLIC :: Filter
  PUBLIC :: Norm
  PUBLIC :: Dot

  !> The constant @f$ \alpha @f$ for adding 2 quadtrees in place.
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_InPlace_Beta
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha
  REAL(SpAMM_KIND) :: SpAMM_Add_BiTree_2_BiTree_RePlace_Beta
  REAL(SpAMM_KIND) :: SpAMM_Add_Identity_2_QuTree_InPlace_Alpha

  !> The SpAMM threshold for the multiply.
  REAL(SpAMM_KIND) :: SpAMM_Threshold_Multiply_QuTree_x_BiTree

  !> Interface for multiplication operations between different SpAMM types.
  !!
  !! The Sparse Approximate Matrix-Multiply (SpAMM):
  !! @f$ C \leftarrow A \cdot B @f$.
  INTERFACE Multiply
     !--------------------------------------------------------
!     MODULE PROCEDURE SpAMM_Multiply_Tree_2d_x_Tree_2d_x_Tree_2d
!     MODULE PROCEDURE SpAMM_Multiply_Tree_2d_x_Tree_2d
!     MODULE PROCEDURE SpAMM_Multiply_Tree_2d_x_Scalar
!     MODULE PROCEDURE SpAMM_Multiply_Tree_2d_x_Tree_1d
!     MODULE PROCEDURE SpAMM_Multiply_Tree_1d_x_Tree_1d
!     MODULE PROCEDURE SpAMM_Multiply_Tree_1d_x_Scalar
     !--------------------------------------------------------
     MODULE PROCEDURE SpAMM_Multiply_QuTree_x_QuTree
     MODULE PROCEDURE SpAMM_Multiply_QuTree_x_Scalar
     MODULE PROCEDURE SpAMM_Multiply_QuTree_x_BiTree
     MODULE PROCEDURE SpAMM_Multiply_BiTree_x_Scalar
     !--------------------------------------------------------
     module procedure spamm_multiply_order_1_x_scalar
     module procedure spamm_multiply_2nd_order_x_2nd_order
     module procedure spamm_multiply_2nd_order_x_scalar
     module procedure spamm_multiply_2nd_order_x_1st_order
  END INTERFACE Multiply

  !> Interface for trace operations.
  INTERFACE Trace
     MODULE PROCEDURE SpAMM_Trace_QuTree
     MODULE PROCEDURE SpAMM_Trace_QuTree_Product
     module procedure spamm_trace_2nd_order
     module procedure spamm_trace_2nd_order_product
  END INTERFACE Trace

  !> Interface for additions operations between different SpAMM types.
  INTERFACE Add
!     MODULE PROCEDURE SpAMM_Add_Tree_2d_2_Tree_2d
!     MODULE PROCEDURE SpAMM_Add_Tree_1d_2_Tree_1d
!     MODULE PROCEDURE SpAMM_Add_Scalar_2_Tree_2d
     !--------------------------------------------------------
     MODULE PROCEDURE SpAMM_Add_QuTree_2_QuTree_InPlace
     MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_InPlace
     MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_RePlace
     MODULE PROCEDURE SpAMM_Add_Identity_2_QuTree_InPlace
     !--------------------------------------------------------
     module procedure spamm_add_identity_to_matrix_2nd_order
     module procedure spamm_add_2nd_order_to_2nd_order
  END INTERFACE Add

  !> Interface for filter operations (thresholding of small matrix elements).
  INTERFACE Filter
     MODULE PROCEDURE SpAMM_Filter_QuTree
  END INTERFACE Filter

  !> Interface for norm operations.
  INTERFACE Norm
     MODULE PROCEDURE SpAMM_Norm_Reduce_BiTree
     MODULE PROCEDURE SpAMM_Norm_Reduce_QuTree
     module procedure spamm_norm_reduce_matrix_2nd_order
  END INTERFACE Norm

  !> Interface for dot product operations.
  INTERFACE Dot
     MODULE PROCEDURE SpAMM_Dot_Product_BiTree
  END INTERFACE Dot

CONTAINS



  INCLUDE 'old_spamm_algebra.F90'

end module spamm_algebra
