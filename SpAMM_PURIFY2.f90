!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    This work is derived under the GPL from the MondoSCF suite of programs,  
!    copyright (2004) the Regents of the University of California,
!    material produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!-------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE SpAMM_QU2BCSR
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE DenMatMethods
  USE MondoLogger
  !       
  USE SpAMM_DERIVED
  USE SpAMM_GLOBALS
  USE SpAMM_MNGMENT
  USE SpAMM_ALGEBRA
  USE SpAMM_CONVERT
  USE SpAMM_PROJECT
  !
  TYPE(BCSR) :: GLOBAL_BCSR
CONTAINS
  SUBROUTINE Copy_GLOBAL_BCSR_2_QuTree(qA)
    TYPE(QuTree), POINTER :: qA
    INTEGER               :: Depth
    REAL(SpAMM_DOUBLE)    :: TInitial, TTotal
    TInitial=SpAMM_IPM_GET_TIME()
    Depth=0
    CALL Copy_GLOBAL_BCSR_2_QuTree_Recur(qA,1,SpAMM_PADDED_MATRIX_DIMENSION,1,SpAMM_PADDED_MATRIX_DIMENSION,Depth) 
    qA%Norm=SQRT(Norm(qA))
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"Copy_GLOBAL_BCSR_2_QuTree",40)
  END SUBROUTINE Copy_GLOBAL_BCSR_2_QuTree

  SUBROUTINE Copy_QuTree_2_GLOBAL_BCSR(qA)
    TYPE(QuTree), POINTER :: qA
    INTEGER               :: Depth,AtA,AtB,NN,P,R
    REAL(SpAMM_DOUBLE)    :: TInitial, TTotal
    TInitial=SpAMM_IPM_GET_TIME()
    !
    P=1 
    R=1 
    GLOBAL_BCSR%RowPt%I(1)=1
    DO AtA=1,NAtoms
       DO AtB=1,NAtoms
          NN=BSiz%I(AtA)*BSiz%I(AtB)
          GLOBAL_BCSR%MTrix%D(R:R+NN-1)=0D0
          GLOBAL_BCSR%ColPt%I(P)=AtB
          GLOBAL_BCSR%BlkPt%I(P)=R
          R=R+NN 
          P=P+1 
          GLOBAL_BCSR%RowPt%I(AtA+1)=P        
       ENDDO
    ENDDO
    GLOBAL_BCSR%NBlks=P-1
    GLOBAL_BCSR%NNon0=R-1
    !
    Depth=0
    CALL Copy_QuTree_2_GLOBAL_BCSR_Recur(qA,1,SpAMM_PADDED_MATRIX_DIMENSION,1,SpAMM_PADDED_MATRIX_DIMENSION,Depth) 
    TTotal=SpAMM_IPM_GET_TIME()-TInitial
    CALL SpAMM_Time_Stamp(TTotal,"Copy_QuTree_2_GLOBAL_BCSR",41)
  END SUBROUTINE Copy_QuTree_2_GLOBAL_BCSR

  RECURSIVE SUBROUTINE Copy_GLOBAL_BCSR_2_QuTree_Recur(qA,LeftRow,RghtRow,LeftCol,RghtCol,Depth) 
    TYPE(QuTree), POINTER :: qA
    INTEGER :: Depth
    INTEGER :: I,LeftRow,RghtRow,LeftCol,RghtCol,HalfRow,HalfCol
    INTEGER :: MA,NA,IAts,IStrt,IStop,J,P,JP,N,M
    INTEGER :: Row,Col,BlkRow,BlkCol
    IF(LeftRow>RghtRow)RETURN
    IF(LeftCol>RghtCol)RETURN    
    IF(.NOT.ASSOCIATED(qA))THEN
       CALL NewQuNode(qA)
    ENDIF
    ! Blocks    
    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
       ! Allocate
       IF(.NOT.ALLOCATED(qA%Blok))THEN
          ALLOCATE(qA%Blok(1:SpAMM_BLOCK_SIZE,1:SpAMM_BLOCK_SIZE))
       END IF
       qA%Blok=Zero
       ! Ugly but it works...  
       DO IAts=1,NAtoms
          MA=BSiz%I(IAts)
          IStrt=GLOBAL_BCSR%RowPt%I(IAts)
          IStop=GLOBAL_BCSR%RowPt%I(IAts+1)-1
          IF(IStrt.NE.0.AND.IStop.NE.0)THEN
             DO JP=IStrt,IStop
                J=GLOBAL_BCSR%ColPt%I(JP)
                P=GLOBAL_BCSR%BlkPt%I(JP)
                NA=BSiz%I(J)
                BlkCol=0
                DO N=1,NA 
                   Col=OffS%I(J)+N-1
                   IF(Col>=LeftCol.AND.Col<=RghtCol)THEN
                      BlkCol=Col-LeftCol+1
                      DO M=1,MA
                         Row=OffS%I(IAts)+M-1
                         IF(Row>=LeftRow.AND.Row<=RghtRow)THEN
                            BlkRow=Row-LeftRow+1
                            qA%Blok(BlkRow,BlkCol)=GLOBAL_BCSR%MTrix%D(P+M+(N-1)*MA-1)
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ELSE
       HalfRow=(RghtRow-LeftRow)/2
       HalfCol=(RghtCol-LeftCol)/2
       CALL Copy_GLOBAL_BCSR_2_QuTree_Recur(qA%Quad00,LeftRow,LeftRow+HalfRow,LeftCol,LeftCol+HalfCol, &
            Depth+1) 
       CALL Copy_GLOBAL_BCSR_2_QuTree_Recur(qA%Quad01,LeftRow,LeftRow+HalfRow,LeftCol+HalfCol+1,RghtCol, &
            Depth+1) 
       CALL Copy_GLOBAL_BCSR_2_QuTree_Recur(qA%Quad10,LeftRow+HalfRow+1,RghtRow,LeftCol,LeftCol+HalfCol, &
            Depth+1) 
       CALL Copy_GLOBAL_BCSR_2_QuTree_Recur(qA%Quad11,LeftRow+HalfRow+1,RghtRow,LeftCol+HalfCol+1,RghtCol, &
            Depth+1) 
    ENDIF
    !
  END SUBROUTINE Copy_GLOBAL_BCSR_2_QuTree_Recur

  RECURSIVE SUBROUTINE Copy_QuTree_2_GLOBAL_BCSR_Recur(qA,LeftRow,RghtRow,LeftCol,RghtCol,Depth) 
    TYPE(QuTree), POINTER :: qA
    INTEGER :: Depth
    INTEGER :: I,LeftRow,RghtRow,LeftCol,RghtCol,HalfRow,HalfCol
    INTEGER :: MA,NA,IAts,IStrt,IStop,J,P,JP,N,M
    INTEGER :: Row,Col,BlkRow,BlkCol
    IF(LeftRow>RghtRow)RETURN
    IF(LeftCol>RghtCol)RETURN    
    IF(.NOT.ASSOCIATED(qA))RETURN
    ! Blocks    
    IF(Depth==SpAMM_TOTAL_DEPTH)THEN
       ! Ugly but it works...  
       DO IAts=1,NAtoms
          MA=BSiz%I(IAts)
          IStrt=GLOBAL_BCSR%RowPt%I(IAts)
          IStop=GLOBAL_BCSR%RowPt%I(IAts+1)-1
          IF(IStrt.NE.0.AND.IStop.NE.0)THEN
             DO JP=IStrt,IStop
                J=GLOBAL_BCSR%ColPt%I(JP)
                P=GLOBAL_BCSR%BlkPt%I(JP)
                NA=BSiz%I(J)
                BlkCol=0
                DO N=1,NA 
                   Col=OffS%I(J)+N-1
                   IF(Col>=LeftCol.AND.Col<=RghtCol)THEN
                      BlkCol=Col-LeftCol+1
                      DO M=1,MA
                         Row=OffS%I(IAts)+M-1
                         IF(Row>=LeftRow.AND.Row<=RghtRow)THEN
                            BlkRow=Row-LeftRow+1
                            GLOBAL_BCSR%MTrix%D(P+M+(N-1)*MA-1)=qA%Blok(BlkRow,BlkCol)
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ELSE
       HalfRow=(RghtRow-LeftRow)/2
       HalfCol=(RghtCol-LeftCol)/2
       CALL Copy_QuTree_2_GLOBAL_BCSR_Recur(qA%Quad00,LeftRow,LeftRow+HalfRow,LeftCol,LeftCol+HalfCol, &
            Depth+1) 
       CALL Copy_QuTree_2_GLOBAL_BCSR_Recur(qA%Quad01,LeftRow,LeftRow+HalfRow,LeftCol+HalfCol+1,RghtCol, &
            Depth+1) 
       CALL Copy_QuTree_2_GLOBAL_BCSR_Recur(qA%Quad10,LeftRow+HalfRow+1,RghtRow,LeftCol,LeftCol+HalfCol, &
            Depth+1) 
       CALL Copy_QuTree_2_GLOBAL_BCSR_Recur(qA%Quad11,LeftRow+HalfRow+1,RghtRow,LeftCol+HalfCol+1,RghtCol, &
            Depth+1) 
    ENDIF
    !
  END SUBROUTINE Copy_QuTree_2_GLOBAL_BCSR_Recur
END MODULE SpAMM_QU2BCSR

PROGRAM DMP_SP2 ! Density matrix purification, SP2 variation
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE DenMatMethods
  USE MondoLogger
  USE SpAMM_QU2BCSR

  IMPLICIT NONE
  !-------------------------------------------------------------------------------
  ! Trace Setting SP2
  !-------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Ne,Lambda, idempotency_error,Occ0,Occ1,Occ2,Occ3
  INTEGER                        :: I, I2, MM, Imin
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='SP2'
  TYPE(BCSR)                     :: F,FT,P,PT,Pold,Tmp1,Tmp2
  TYPE(DBL_RNK2)                 :: dF
  TYPE(DBL_VECT)                 :: EigenV
  !-------------------------------------------------------------------------------
  TYPE(QuTree),POINTER           :: qF=>NULL(), qP=>NULL(), qTmp1=>NULL(), qTmp2=>NULL()
  REAL(SpAMM_KIND)               :: TrE,HalfNe,Nocc
  !-------------------------------------------------------------------------------
#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  !
  CALL New(GLOBAL_BCSR,N_O=(/1+NAtoms,1+NAtoms**2,1+NBasF**2/))
  !
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
    CALL Get(GLOBAL_BCSR,FFile)
  ELSE
    CALL Get(GLOBAL_BCSR,TrixFile('OrthoF',Args,0))
  ENDIF

  ! Following call raises all kinds of questions about spectral 
  ! remaping, b/c using perfect idempotent matrix still leads to 
  ! lots of purification work. Shouldnt be ...
  !
!   CALL Get(GLOBAL_BCSR,TrixFile('OrthoD',Args,0))
!!$
!!$  CALL New(dF,(/NBasF,NBasF/))
!!$  CALL SetEq(dF,GLOBAL_BCSR)
!!$  CALL New(EigenV,NBasF)
!!$  CALL SetEq(EigenV,Zero)
!!$  CALL MDiag_DSYEVD(dF,NBasF,EigenV,0)
!!$
!!$  WRITE(*,*)'EXTREME EIGENVALUES = ',EigenV%D

!  STOP
!  CALL Print_CheckSum_BCSR(GLOBAL_BCSR,' OrthoF(A) in ',Unit_O=6) 




  !-----------------------------------------------------------
  CALL SpAMM_Init_Globals(NBasF,1)
  CALL NewQuNode(qF,init=.TRUE.)
  CALL Copy_GLOBAL_BCSR_2_QuTree(qF)   
  CALL Copy(qF,qP)
  CALL NewQuNode(qTmp1,init=.TRUE.)
  CALL AllocateFull(qTmp1)
  !-----------------------------------------------------------
  CALL RemapSpectralBounds201(qP)
!!$
!!$  CALL Copy_QuTree_2_GLOBAL_BCSR(qP) 
!!$  CALL New(dF,(/NBasF,NBasF/))
!!$  CALL SetEq(dF,GLOBAL_BCSR)
!!$  CALL New(EigenV,NBasF)
!!$  CALL SetEq(EigenV,Zero)
!!$  CALL MDiag_DSYEVD(dF,NBasF,EigenV,0)
!!$
!!$  WRITE(*,*)'EXTREME EIGENVALUES = ',EigenV%D

!  STOP


  HalfNe=SpAMM_Half*FLOAT(NEl)
  Occ0=Zero
  Occ1=Zero
  Occ2=Zero
  Occ3=Zero
  Imin = 20
  DO I=1,100
     CALL SpAMM_TC2(qP,qTmp1,HalfNe,Nocc)
     WRITE(*,*)' I ',I,ABS(Nocc*SpAMM_Two),ABS(Nocc*SpAMM_Two-FLOAT(NEl))/FLOAT(NEl)

     Occ0=Nocc
     IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) THEN
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "converged in "//TRIM(IntToChar(I))//" iterations")
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "Idempotency error = "//TRIM(DblToChar(ABS(Occ0-Occ1))))
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "Previous idempotency error = "//TRIM(DblToChar(ABS(Occ2-Occ3))))
        EXIT
     ENDIF
     Occ3 = Occ2
     Occ2 = Occ1
     Occ1 = Occ0
  ENDDO
  !-----------------------------------------------------------
  CALL Multiply(qP,qF,qTmp1)
  TrE=Trace(qTmp1)
  CALL Copy_QuTree_2_GLOBAL_BCSR(qP) 
  CALL SpAMM_Time_Stamp()
  WRITE(*,*)' E_el = ',TrE
  !-----------------------------------------------------------
  ! Orthogonal put and xform to AO rep and put
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL PutXForm(Prog,Args,GLOBAL_BCSR,Tmp1,Tmp2)
  !
  CALL Delete(GLOBAL_BCSR)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  !
  CALL ShutDown(Prog)

END PROGRAM DMP_SP2
