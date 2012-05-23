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
!    PACKAGE FOR THE SPARSE APPROXIMATE MATRIX MULTIPLY (SpAMMpack)
!    Matt Challacombe and Nick Bock
!------------------------------------------------------------------------------

!> @brief
!! Defines operations on SpAMM trees.
MODULE SpAMM_ALGEBRA
  USE  SpAMM_DERIVED
  USE  SpAMM_GLOBALS
  USE  SpAMM_MNGMENT

  IMPLICIT NONE

  !===============================================================================
  !  GLOBALS (we don't want these on the stack ...)
  !===============================================================================

  !> The constant @f$ \alpha @f$ for adding 2 quadtrees in place.
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_QuTree_2_QuTree_InPlace_Alpha
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_QuTree_2_QuTree_InPlace_Beta
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_BiTree_2_BiTree_InPlace_Beta
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_BiTree_2_BiTree_RePlace_Beta
  REAL(SpAMM_DOUBLE) :: SpAMM_Add_Identity_2_QuTree_InPlace_Alpha
  REAL(SpAMM_DOUBLE) :: SpAMM_Threshold_Multiply_QuTree_x_QuTree
  REAL(SpAMM_DOUBLE) :: SpAMM_Threshold_Multiply_QuTree_x_BiTree

  !===============================================================================
  !  EXTERNAL
  !===============================================================================
#ifdef SPAMM_DOUBLE
  EXTERNAL :: DGEMM
#else
  EXTERNAL :: SGEMM
#endif

  !===============================================================================
  !  INTERFACE BLOCKS
  !===============================================================================

  !> @brief
  !! Interface for multiplication operations between different SpAMM types.
  !!
  !! @details
  !! The Sparse Approximate Matrix-Multiply (SpAMM):
  !! @f$ C \leftarrow A \times B @f$.
  INTERFACE Multiply
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_QuTree
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_Scalar
    MODULE PROCEDURE SpAMM_Multiply_QuTree_x_BiTree
    MODULE PROCEDURE SpAMM_Multiply_BiTree_x_Scalar
  END INTERFACE

  !> @brief
  !! Interface for trace operations.
  INTERFACE Trace
    MODULE PROCEDURE SpAMM_Trace_QuTree
    MODULE PROCEDURE SpAMM_Trace_QuTree_Product
  END INTERFACE

  !> @brief
  !! Interface for additions operations between different SpAMM types.
  INTERFACE Add
    MODULE PROCEDURE SpAMM_Add_QuTree_2_QuTree_InPlace
    MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_InPlace
    MODULE PROCEDURE SpAMM_Add_BiTree_2_BiTree_RePlace
    MODULE PROCEDURE SpAMM_Add_Identity_2_QuTree_InPlace
  END INTERFACE

  !> @brief
  !! Interface for filter operations (thresholding of small matrix elements).
  INTERFACE Filter
    MODULE PROCEDURE SpAMM_Filter_QuTree
  END INTERFACE

  !> @brief
  !! Interface for norm operations.
  INTERFACE Norm
    MODULE PROCEDURE SpAMM_Norm_Reduce_BiTree
    MODULE PROCEDURE SpAMM_Norm_Reduce_QuTree
  END INTERFACE

  !> @brief
  !! Interface for dot product operations.
  INTERFACE Dot
    MODULE PROCEDURE SpAMM_Dot_Product_BiTree
  END INTERFACE

CONTAINS

  !> Multiplication operation between a quadtree and a quadtree.
  SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree(qA,qB,qC,LocalThreshold)

    TYPE(QuTree), POINTER              :: qA,qB
    TYPE(QuTree), POINTER, INTENT(OUT) :: qC
    INTEGER :: Depth
    REAL(SpAMM_KIND), OPTIONAL         :: LocalThreshold
    REAL(SpAMM_DOUBLE)                 :: TInitial, TTotal

    IF(PRESENT(LocalThreshold))THEN
      SpAMM_Threshold_Multiply_QuTree_x_QuTree=LocalThreshold
    ELSE
      SpAMM_Threshold_Multiply_QuTree_x_QuTree=SpAMM_PRODUCT_TOLERANCE
    ENDIF
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    CALL SpAMM_Multiply_QuTree_x_Scalar(qC,SpAMM_Zero)
    !$OMP TASK UNTIED SHARED(qA,qB,qC)
    CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC,qA,qB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_QuTree",1)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! The Sparse Approximate Matrix-Multiply (SpAMM): D <- A.B.C
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Scalar multiply: A <- a*A
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar(qA,a)

    INTEGER :: Depth
    TYPE(QuTree), POINTER    :: qA
    REAL(SpAMM_KIND)         :: a
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    IF(.NOT.ASSOCIATED(qA))RETURN
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK SHARED(qA)
    CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA,a,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_Scalar",3)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree In Place Add: A <- a*A + b*B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace(qA,qB,Alpha,Beta)

    TYPE(QuTree), POINTER, INTENT(INOUT) :: qA
    TYPE(QuTree), POINTER, INTENT(IN)    :: qB
    REAL(SpAMM_KIND), OPTIONAL :: Alpha,Beta
    INTEGER                    :: Depth
    REAL(SpAMM_KIND)           :: Norm
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    IF(PRESENT(Alpha))THEN
      SpAMM_Add_QuTree_2_QuTree_InPlace_Alpha=Alpha
    ELSE
      SpAMM_Add_QuTree_2_QuTree_InPlace_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
      SpAMM_Add_QuTree_2_QuTree_InPlace_Beta=Beta
    ELSE
      SpAMM_Add_QuTree_2_QuTree_InPlace_Beta=SpAMM_One
    ENDIF
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(qA,qB)
    CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA,qB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_QuTree_2_QuTree_InPlace",4)

  END SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree In Place Add: A <- A + a*I
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace(qA,Alpha)

    TYPE(QuTree), POINTER                :: qA
    REAL(SpAMM_KIND)                     :: Alpha
    INTEGER                              :: Depth
    REAL(SpAMM_KIND)                     :: Norm
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    SpAMM_Add_Identity_2_QuTree_InPlace_Alpha=Alpha
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA,1,SpAMM_PADDED_MATRIX_DIMENSION,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_QuTree_2_QuTree_InPlace",5)

  END SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Trace for QuTree: a = trace[A]
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FUNCTION SpAMM_Trace_QuTree(qA) RESULT(a)

    TYPE(QuTree), POINTER  :: qA
    REAL(SpAMM_KIND)       :: a
    INTEGER                :: Depth
    REAL(SpAMM_DOUBLE)     :: TInitial, TTotal
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(qA,a)
    a=SpAMM_Trace_QuTree_Recur(qA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree",6)

  END FUNCTION SpAMM_Trace_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Trace for QuTree: a = trace[A.B]
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FUNCTION SpAMM_Trace_QuTree_Product(qA,qB) RESULT(a)

    TYPE(QuTree), POINTER  :: qA,qB
    REAL(SpAMM_KIND)       :: a
    INTEGER                :: Depth
    REAL(SpAMM_DOUBLE)     :: TInitial, TTotal
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(qA,a)
    a=SpAMM_Trace_QuTree_Product_Recur(qA,qB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Trace_QuTree",7)

  END FUNCTION SpAMM_Trace_QuTree_Product

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Filter for the QuTree: \tilde{A}=filter[A,tau]
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Filter_QuTree(qA,Tau)

    TYPE(QuTree), POINTER  :: qA
    REAL(SpAMM_KIND)       :: Tau
    INTEGER                :: Depth
    REAL(SpAMM_DOUBLE)     :: TInitial, TTotal
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(qA)
    CALL SpAMM_Filter_QuTree_Recur(qA,Tau,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Filter_QuTree",8)

  END SUBROUTINE SpAMM_Filter_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! L_2 norm for QuTrees
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FUNCTION SpAMM_Norm_Reduce_QuTree(qA) RESULT(Norm)

    INTEGER :: Depth
    TYPE(QuTree),POINTER :: qA
    REAL(SpAMM_KIND) :: Norm
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK SHARED(Norm,qA)
    Norm=SpAMM_Norm_Reduce_QuTree_Recur(qA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Norm_Reduce_QuTree",9)

  END FUNCTION SpAMM_Norm_Reduce_QuTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! C <- Alpha*A+Beta*B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace(bC,Alpha,bA,Beta,bB)

    TYPE(BiTree), POINTER, INTENT(IN)    :: bA,bB
    TYPE(BiTree), POINTER, INTENT(INOUT) :: bC
    REAL(SpAMM_KIND), OPTIONAL           :: Alpha,Beta
    INTEGER                              :: Depth
    REAL(SpAMM_DOUBLE)                   :: TInitial, TTotal
    Depth=0
    IF(PRESENT(Alpha))THEN
      SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha=Alpha
    ELSE
      SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
      SpAMM_Add_BiTree_2_BiTree_RePlace_Beta=Beta
    ELSE
      SpAMM_Add_BiTree_2_BiTree_RePlace_Beta=SpAMM_One
    ENDIF
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(bA,bB)
    CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA,bB,bC,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_BiTree_2_BiTree_RePlace",10)

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! A<-Apha*A+Beta*B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace(bA,Alpha,bB,Beta)

    TYPE(BiTree), POINTER, INTENT(INOUT) :: bA
    TYPE(BiTree), POINTER, INTENT(IN)    :: bB
    REAL(SpAMM_KIND), OPTIONAL           :: Alpha,Beta
    INTEGER                              :: Depth
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    IF(PRESENT(Alpha))THEN
      SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha=Alpha
    ELSE
      SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha=SpAMM_One
    ENDIF
    IF(PRESENT(Beta))THEN
      SpAMM_Add_BiTree_2_BiTree_InPlace_Beta=Beta
    ELSE
      SpAMM_Add_BiTree_2_BiTree_InPlace_Beta=SpAMM_One
    ENDIF
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK UNTIED SHARED(bA,bB)
    CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA,bB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Add_BiTree_2_BiTree_InPlace",11)

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Inner product: (A,B)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FUNCTION SpAMM_Dot_Product_BiTree(bA,bB) RESULT(dot)

    INTEGER              :: Depth
    TYPE(BiTree),POINTER :: bA,bB
    REAL(SpAMM_KIND)     :: dot
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK SHARED(dot,bA,bB)
    dot=SpAMM_Dot_Product_BiTree_Recur(bA,bB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Dot_Product_BiTree",12)

  END FUNCTION SpAMM_Dot_Product_BiTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Sparse Approximate Matrix-Vector Multiply (SpAMV): C <- A.B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree(qA,bB,bC,LocalThreshold)

    TYPE(QuTree), POINTER     :: qA
    TYPE(BiTree), POINTER     :: bB,bC
    INTEGER                   :: Depth
    REAL(SpAMM_KIND),OPTIONAL :: LocalThreshold
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    IF(PRESENT(LocalThreshold))THEN
      SpAMM_Threshold_Multiply_QuTree_x_BiTree=LocalThreshold
    ELSE
      SpAMM_Threshold_Multiply_QuTree_x_BiTree=SpAMM_PRODUCT_TOLERANCE
    ENDIF
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    CALL SpAMM_Multiply_BiTree_x_Scalar(bC,SpAMM_Zero)
    !$OMP TASK UNTIED SHARED(qA,bB,bC)
    CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC,qA,bB,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_QuTree_x_BiTree",13)

  END SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Scalar multiply: A <- a*A
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar(bA,a)

    INTEGER :: Depth
    TYPE(BiTree), POINTER    :: bA
    REAL(SpAMM_KIND)         :: a
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    IF(.NOT.ASSOCIATED(bA))RETURN
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK SHARED(bA)
    CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA,a,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Multiply_BiTree_x_Scalar",14)

  END SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Norm for BiTrees
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FUNCTION SpAMM_Norm_Reduce_BiTree(bA) RESULT(Norm)

    INTEGER              :: Depth
    TYPE(BiTree),POINTER :: bA
    REAL(SpAMM_KIND)     :: Norm
    REAL(SpAMM_DOUBLE)                                  :: TInitial, TTotal
    Depth=0
    TInitial=SpAMM_IPM_GET_TIME()
    !$OMP TASK SHARED(Norm,bA)
    Norm=SpAMM_Norm_Reduce_BiTree_Recur(bA,Depth)
    !$OMP END TASK
    !$OMP TASKWAIT
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"SpAMM_Norm_Reduce_BiTree",15)

  END FUNCTION SpAMM_Norm_Reduce_BiTree

  !=================================================================
  ! RECURSIVE LINEAR ALGEBRA ROUTINES
  !=================================================================
  ! The recursive Sparse Approximate Matrix-Multiply (SpAMM)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur(qC,qA,qB,Depth)

    TYPE(QuTree), POINTER :: qC,qA,qB
    INTEGER :: Depth
    LOGICAL :: DepthOK,Go_00x00,Go_00x01,Go_10x00,Go_10x01,Go_01x10,Go_01x11,Go_11x10,Go_11x11
    ! Associated
    IF(ASSOCIATED(qA).AND.ASSOCIATED(qB))THEN
      ! Estimate
      IF(qA%Norm*qB%Norm<SpAMM_Threshold_Multiply_QuTree_x_QuTree)RETURN
      IF(.NOT.ASSOCIATED(qC))THEN
        !$OMP CRITICAL
        ALLOCATE(qC)
        !$OMP END CRITICAL
      ENDIF
      ! Blocks
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        ! Allocate
        IF(.NOT.ALLOCATED(qC%Blok))THEN
          !$OMP CRITICAL
          ALLOCATE(qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
          qC%Blok=SpAMM_Zero
          !$OMP END CRITICAL
        END IF
        ! Accumulate
        qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE)=         &
          qC%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE)+MATMUL(  &
          qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE),         &
          qB%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
      ELSE
#ifdef _OPENMP
        ! Put a check on the stack
        DepthOK=.TRUE.
        !          DepthOK=MOD(Depth,2)==0
        Go_00x00=DepthOK.AND.qA%Quad00%Norm*qB%Quad00%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_00x01=DepthOK.AND.qA%Quad00%Norm*qB%Quad01%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_10x00=DepthOK.AND.qA%Quad10%Norm*qB%Quad00%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_10x01=DepthOK.AND.qA%Quad10%Norm*qB%Quad01%Norm>SpAMM_RECURSION_NORMD_CUTOFF
#endif
        ! 00=00*00
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_00x00)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad00,qA%Quad00,qB%Quad00,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_00x01)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad01,qA%Quad00,qB%Quad01,Depth+1)
        !$OMP END TASK
        ! 10=10*00
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_10x00)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad10,qA%Quad10,qB%Quad00,Depth+1)
        !$OMP END TASK
        ! 11=10*01
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_10x01)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad11,qA%Quad10,qB%Quad01,Depth+1)
        !$OMP END TASK
#ifdef _OPENMP
        Go_01x10=DepthOK.AND.qA%Quad01%Norm*qB%Quad10%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_01x11=DepthOK.AND.qA%Quad01%Norm*qB%Quad11%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_11x10=DepthOK.AND.qA%Quad11%Norm*qB%Quad10%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_11x11=DepthOK.AND.qA%Quad11%Norm*qB%Quad11%Norm>SpAMM_RECURSION_NORMD_CUTOFF
#endif
        !$OMP TASKWAIT
        ! 00=00*00+01*10
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_01x10)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad00,qA%Quad01,qB%Quad10,Depth+1)
        !$OMP END TASK
        ! 01=00*01+01*11
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_01x11)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad01,qA%Quad01,qB%Quad11,Depth+1)
        !$OMP END TASK
        ! 10=10*00+11*10
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_11x10)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad10,qA%Quad11,qB%Quad10,Depth+1)
        !$OMP END TASK
        ! 11=10*01+11*11
        !$OMP TASK UNTIED SHARED(qA,qB,qC) IF(Go_11x11)
        CALL SpAMM_Multiply_QuTree_x_QuTree_Recur(qC%Quad11,qA%Quad11,qB%Quad11,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
      ENDIF
    ENDIF
  END SUBROUTINE SpAMM_Multiply_QuTree_x_QuTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur(qA,a,Depth)

    INTEGER :: Depth
    TYPE(QuTree), POINTER   :: qA
    REAL(SpAMM_KIND)        :: a
    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(qA%Blok))THEN
      qA%Norm=qA%Norm*ABS(a)
      qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE)=qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE)*a
    ELSE
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad00,a,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED  SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad01,a,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad10,a,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Multiply_QuTree_x_Scalar_Recur(qA%Quad11,a,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      qA%Norm=qA%Norm*ABS(a)
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_Scalar_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree In Place Add: A <- alpha*A + beta*B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA,qB,Depth)

    TYPE(QuTree),POINTER :: qA,qB
    INTEGER              :: Depth
    LOGICAL              :: TA, TB
    TA=ASSOCIATED(qA)
    TB=ASSOCIATED(qB)
    IF(TA.AND.TB)THEN
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        qA%Blok=SpAMM_Add_QuTree_2_QuTree_InPlace_Alpha*qA%Blok &
          +SpAMM_Add_QuTree_2_QuTree_InPlace_Beta *qB%Blok
      ELSE
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad00,qB%Quad00,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad01,qB%Quad01,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad10,qB%Quad10,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_QuTree_2_QuTree_InPlace_Recur(qA%Quad11,qB%Quad11,Depth+1)
        !$OMP END TASK
      ENDIF
      !$OMP TASKWAIT !! << WTF IS THIS TASKWAIT IMPORTANT FOR IFORT??
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(qA,qB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_QuTree_2_QuTree_Recur(qB,qA,Depth)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Add_QuTree_2_QuTree_InPlace_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree In Place Add: A <- A + alpha*I
  ! Note, this routine is not empowered to deal with
  ! case of missing diagonal blocks from the in place QuTree.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA,Left,Rght,Depth)

    TYPE(QuTree),POINTER :: qA
    INTEGER              :: Depth
    INTEGER              :: Left,Rght,Half,I
    real*8 ss
    IF(Left>SpAMM_MATRIX_DIMENSION.OR.Left>Rght)THEN
      RETURN
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
      IF(.NOT.ASSOCIATED(qA))THEN
        CALL NewQuNode(qA)
        ALLOCATE(qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
        qA%Blok=SpAMM_Zero
      ENDIF
      DO I=1,MIN(SpAMM_BLOCK_SIZE,SpAMM_MATRIX_DIMENSION-Left+1)
        qA%Blok(I,I)=qA%Blok(I,I)+SpAMM_Add_Identity_2_QuTree_InPlace_Alpha
      ENDDO
    ELSE
      Half=(Rght-Left)/2
      !$OMP TASK UNTIED SHARED(qA) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA%Quad00,Left,Left+Half,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Add_Identity_2_QuTree_InPlace_Recur(qA%Quad11,Left+Half+1,Rght,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT !! << WTF IS THIS TASKWAIT IMPORTANT FOR IFORT??
    ENDIF

  END SUBROUTINE SpAMM_Add_Identity_2_QuTree_InPlace_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree recursive trace
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE FUNCTION SpAMM_Trace_QuTree_Recur(qA,Depth) RESULT(Trace)

    TYPE(QuTree), POINTER  :: qA
    REAL(SpAMM_KIND) :: Trace,Trace00,Trace11
    INTEGER :: Depth, I
    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      Trace=SpAMM_Zero
      IF(.NOT.ASSOCIATED(qA))RETURN
      DO I=1,SpAMM_BLOCK_SIZE
        Trace=Trace+qA%Blok(I,I)
      ENDDO
    ELSEIF(.NOT.ASSOCIATED(qA%Quad00).AND. &
        .NOT.ASSOCIATED(qA%Quad11))THEN
      Trace=SpAMM_Zero
    ELSEIF(.NOT.ASSOCIATED(qA%Quad11))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace=SpAMM_Trace_QuTree_Recur(qA%Quad00,Depth+1)
      !$OMP END TASK
    ELSEIF(.NOT.ASSOCIATED(qA%Quad00))THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace=SpAMM_Trace_QuTree_Recur(qA%Quad11,Depth+1)
      !$OMP END TASK
    ELSE
      !$OMP TASK UNTIED SHARED(qA,Trace00) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace00=SpAMM_Trace_QuTree_Recur(qA%Quad00,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Trace11) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Trace11=SpAMM_Trace_QuTree_Recur(qA%Quad11,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      Trace=Trace00+Trace11
    ENDIF

  END FUNCTION SpAMM_Trace_QuTree_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree recursive trace
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE FUNCTION SpAMM_Trace_QuTree_Product_Recur(qA,qB,Depth) RESULT(Trace)

    TYPE(QuTree), POINTER  :: qA,qB
    REAL(SpAMM_KIND) :: Trace,Trace00,Trace11
    INTEGER :: Depth, I,J
    Trace=SpAMM_Zero
    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(.NOT.ASSOCIATED(qB))RETURN
    IF(qA%Norm*qB%Norm<SpAMM_Threshold_Multiply_QuTree_x_QuTree)RETURN
    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
      DO I=1,SpAMM_BLOCK_SIZE
        Trace=Trace+DOT_PRODUCT(qA%Blok(I,1:SpAMM_BLOCK_SIZE),qB%Blok(1:SpAMM_BLOCK_SIZE,I))
      ENDDO
    ELSE
      ! 00=00*00
      !$OMP TASK UNTIED SHARED(qA,qB,Trace00)
      Trace00=SpAMM_Trace_QuTree_Product_Recur(qA%Quad00,qB%Quad00,Depth+1)
      !$OMP END TASK
      ! 11=10*01
      !$OMP TASK UNTIED SHARED(qA,qB,Trace11)
      Trace11=SpAMM_Trace_QuTree_Product_Recur(qA%Quad10,qB%Quad01,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      ! 00=00*00+01*10
      !$OMP TASK UNTIED SHARED(qA,qB,Trace00)
      Trace00=Trace00+SpAMM_Trace_QuTree_Product_Recur(qA%Quad01,qB%Quad10,Depth+1)
      !$OMP END TASK
      ! 11=10*01+11*11
      !$OMP TASK UNTIED SHARED(qA,qB,Trace11)
      Trace11=Trace11+SpAMM_Trace_QuTree_Product_Recur(qA%Quad11,qB%Quad11,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      Trace=Trace00+Trace11
    ENDIF

  END FUNCTION SpAMM_Trace_QuTree_Product_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! QuTree filter: \tilde{A}=filter[A,tau]
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Filter_QuTree_Recur(qA,Tau,Depth)

    TYPE(QuTree), POINTER  :: qA
    REAL(SpAMM_KIND)       :: Tau
    INTEGER                :: Depth
    IF(.NOT.ASSOCIATED(qA))RETURN
    IF(qA%Norm<Tau)THEN
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Delete_QuTree_Recur(qA,Depth)
      !$OMP END TASK
      !$OMP TASKWAIT
      DEALLOCATE(qA)
    ELSE
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad00,Tau,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad01,Tau,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad10,Tau,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Filter_QuTree_Recur(qA%Quad11,Tau,Depth+1)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Filter_QuTree_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! L_2 norm for QuTrees
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE FUNCTION SpAMM_Norm_Reduce_QuTree_Recur(qA,Depth) RESULT(Norm)

    TYPE(QuTree), POINTER :: qA
    INTEGER               :: Depth
    REAL(SpAMM_KIND)      :: Norm,Norm00,Norm01,Norm10,Norm11
    IF(.NOT.ASSOCIATED(qA))THEN
      Norm=SpAMM_Zero
      RETURN
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
      !qA%Siz==SpAMM_BLOCK_SIZE.AND.ALLOCATED(qA%Blok))THEN
      Norm=SUM(qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE)**2)
      qA%Norm=SQRT(Norm)
    ELSE
      !$OMP TASK UNTIED SHARED(qA,Norm00) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm00=SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad00,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Norm01) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm01=SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad01,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Norm10) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm10=SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad10,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(qA,Norm11) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm11=SpAMM_Norm_Reduce_QuTree_Recur(qA%Quad11,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      Norm=Norm00+Norm01+Norm10+Norm11
      qA%Norm=SQRT(Norm)
    ENDIF

  END FUNCTION SpAMM_Norm_Reduce_QuTree_Recur

  !=================================================================
  ! RECURSIVE LINEAR ALGEBRA ROUTINES ON ROW TREE VECTORS
  !=================================================================
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree_Recur(bC,qA,bB,Depth)

    TYPE(QuTree), POINTER :: qA
    TYPE(BiTree), POINTER :: bB,bC
    INTEGER               :: Depth
    LOGICAL               :: DepthOK,Go_00x0,Go_01x1,Go_10x0,Go_11x1
    ! Associated
    IF(ASSOCIATED(qA).AND.ASSOCIATED(bB))THEN
      ! Estimate
      IF(qA%Norm*bB%Norm<SpAMM_Threshold_Multiply_QuTree_x_BiTree)RETURN
      IF(.NOT.ASSOCIATED(bC))THEN
        !$OMP CRITICAL
        ALLOCATE(bC)
        !$OMP END CRITICAL
      ENDIF
      ! Blocks
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        ! Allocate
        IF(.NOT.ALLOCATED(bC%Vect))THEN
          !$OMP CRITICAL
          ALLOCATE(bC%Vect(1:SpAMM_BLOCK_SIZE))
          bC%Vect=SpAMM_Zero
          !$OMP END CRITICAL
        END IF
        ! Accumulate
        bC%Vect(1:SpAMM_BLOCK_SIZE)=bC%Vect(1:SpAMM_BLOCK_SIZE)+MATMUL( &
          qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE),bB%Vect(1:SpAMM_BLOCK_SIZE))
      ELSE
#ifdef _OPENMP
        ! Put a check on the stack
        DepthOK=.TRUE.
        ! DepthOK=MOD(Depth,2)==0
        Go_00x0=DepthOK.AND.qA%Quad00%Norm*bB%Sect0%Norm>SpAMM_RECURSION_NORMD_CUTOFF
          Go_01x1=DepthOK.AND.qA%Quad01%Norm*bB%Sect1%Norm>SpAMM_RECURSION_NORMD_CUTOFF
        Go_10x0=DepthOK.AND.qA%Quad10%Norm*bB%Sect0%Norm>SpAMM_RECURSION_NORMD_CUTOFF
          Go_11x1=DepthOK.AND.qA%Quad11%Norm*bB%Sect1%Norm>SpAMM_RECURSION_NORMD_CUTOFF
#endif
        ! 0=00*0
        !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(Go_00x0)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect0,qA%Quad00,bB%Sect0,Depth+1)
        !$OMP END TASK
        ! 1=10*0
        !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(Go_10x0)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect1,qA%Quad10,bB%Sect0,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
        ! 0=00*0+01*1
        !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(Go_01x1)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect0,qA%Quad01,bB%Sect1,Depth+1)
        !$OMP END TASK
        ! 1=10*0+11*1
          !$OMP TASK UNTIED SHARED(qA,bB,bC) IF(Go_11x1)
        CALL SpAMM_Multiply_QuTree_x_BiTree_Recur(bC%Sect1,qA%Quad11,bB%Sect1,Depth+1)
        !$OMP END TASK
        !$OMP TASKWAIT
      ENDIF
    ENDIF

  END SUBROUTINE SpAMM_Multiply_QuTree_x_BiTree_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! BiTree in place add: A <- alpha*A + beta*B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA,bB,Depth)

    TYPE(BiTree),POINTER :: bA,bB
    INTEGER              :: Depth
    LOGICAL              :: TA, TB
    TA=ASSOCIATED(bA)
    TB=ASSOCIATED(bB)
    IF(TA.AND.TB)THEN
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        bA%Vect=SpAMM_Add_BiTree_2_BiTree_InPlace_Alpha*bA%Vect &
          +SpAMM_Add_BiTree_2_BiTree_InPlace_Beta *bB%Vect
      ELSE
        !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA%Sect0,bB%Sect0,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_InPlace_Recur(bA%Sect1,bB%Sect1,Depth+1)
        !$OMP END TASK
      ENDIF
      !$OMP TASKWAIT !! << WTF IS THIS TASKWAIT IMPORTANT FOR IFORT??
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(bA,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bB,bA,Depth)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_InPlace_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! BiTree add with replacement: C <- alpha*A + beta*B
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA,bB,bC,Depth)

    TYPE(BiTree),POINTER :: bA,bB,bC
    INTEGER              :: Depth
    LOGICAL              :: TA, TB
    TA=ASSOCIATED(bA)
    TB=ASSOCIATED(bB)
    IF(TA.AND.TB)THEN
      IF(Depth==SpAMM_TOTAL_DEPTH)THEN
        bC%Vect=SpAMM_Add_BiTree_2_BiTree_RePlace_Alpha*bA%Vect &
          +SpAMM_Add_BiTree_2_BiTree_RePlace_Beta *bB%Vect
      ELSE
        !$OMP TASK UNTIED SHARED(bA,bB,bC) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA%Sect0,bB%Sect0,bC%Sect0,Depth+1)
        !$OMP END TASK
        !$OMP TASK UNTIED SHARED(bA,bB,bC) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
        CALL SpAMM_Add_BiTree_2_BiTree_RePlace_Recur(bA%Sect1,bB%Sect1,bC%Sect1,Depth+1)
        !$OMP END TASK
      ENDIF
      !$OMP TASKWAIT !! << WTF IS THIS TASKWAIT IMPORTANT FOR IFORT??
    ELSEIF(.NOT.TA.AND.TB)THEN
      !$OMP TASK UNTIED SHARED(bC,bB) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bB,bC,Depth)
      !$OMP END TASK
    ELSEIF(.NOT.TB.AND.TA)THEN
      !$OMP TASK UNTIED SHARED(bC,bA) IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Copy_BiTree_2_BiTree_Recur(bA,bC,Depth)
      !$OMP END TASK
    ENDIF

  END SUBROUTINE SpAMM_Add_BiTree_2_BiTree_RePlace_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! L_2 norm for BiTrees
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE FUNCTION SpAMM_Norm_Reduce_BiTree_Recur(bA,Depth) RESULT(Norm)

    TYPE(BiTree), POINTER :: bA
    INTEGER               :: Depth
    REAL(SpAMM_KIND)      :: Norm,Norm0,Norm1
    IF(.NOT.ASSOCIATED(bA))THEN
      Norm=SpAMM_Zero
      RETURN
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
      Norm=SUM(bA%Vect(1:SpAMM_BLOCK_SIZE)**2)
      bA%Norm=SQRT(Norm)
    ELSE
       !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm0=SpAMM_Norm_Reduce_BiTree_Recur(bA%Sect0,Depth+1)
      !$OMP END TASK
       !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Norm1=SpAMM_Norm_Reduce_BiTree_Recur(bA%Sect1,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      Norm=Norm0+Norm1
      bA%Norm=SQRT(Norm)
    ENDIF

  END FUNCTION SpAMM_Norm_Reduce_BiTree_Recur

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Dot product for BiTrees
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RECURSIVE FUNCTION SpAMM_Dot_Product_BiTree_Recur(bA,bB,Depth) RESULT(Dot)

    TYPE(BiTree), POINTER :: bA,bB
    INTEGER               :: Depth
    REAL(SpAMM_KIND)      :: Dot, Dot0, Dot1
    IF(.NOT.ASSOCIATED(bA))THEN
      Dot=SpAMM_Zero
      RETURN
    ELSEIF(.NOT.ASSOCIATED(bB))THEN
      Dot=SpAMM_Zero
      RETURN
    ELSEIF(Depth==SpAMM_TOTAL_DEPTH)THEN
      Dot=DOT_PRODUCT(bA%Vect(1:SpAMM_BLOCK_SIZE),bB%Vect(1:SpAMM_BLOCK_SIZE))
    ELSE
      !$OMP TASK UNTIED SHARED(bA,bB,Dot0) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Dot0=SpAMM_Dot_Product_BiTree_Recur(bA%Sect0,bB%Sect0,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED SHARED(bA,bB,Dot1) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      Dot1=SpAMM_Dot_Product_BiTree_Recur(bA%Sect1,bB%Sect1,Depth+1)
      !$OMP END TASK
      !$OMP TASKWAIT
      Dot=Dot0+Dot1
    ENDIF

  END FUNCTION SpAMM_Dot_Product_BiTree_Recur

  RECURSIVE SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar_Recur(bA,a,Depth)

    INTEGER :: Depth
    TYPE(BiTree), POINTER   :: bA
    REAL(SpAMM_KIND)        :: a
    IF(.NOT.ASSOCIATED(bA))RETURN

    IF(Depth==SpAMM_TOTAL_DEPTH.AND.ALLOCATED(bA%Vect))THEN
      bA%Norm=bA%Norm*ABS(a)
      bA%Vect=bA%Vect*a
    ELSE
      !$OMP TASK UNTIED SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA%Sect0,a,Depth+1)
      !$OMP END TASK
      !$OMP TASK UNTIED  SHARED(bA) &
      !$OMP&     IF(Depth<SpAMM_RECURSION_DEPTH_CUTOFF)
      CALL SpAMM_Multiply_BiTree_x_Scalar_Recur(bA%Sect1,a,Depth+1)
      !$OMP END TASK
      bA%Norm=bA%Norm*ABS(a)
    ENDIF

  END SUBROUTINE SpAMM_Multiply_BiTree_x_Scalar_Recur

END MODULE SpAMM_ALGEBRA
