#include "config_fortran.h"

MODULE SpAMMPACK_MANAGEMENT

  USE, INTRINSIC :: iso_C_binding
  USE SpAMMPACK_TYPES

  IMPLICIT NONE

  INTERFACE SpAMM_New

    MODULE PROCEDURE &
        SpAMM_New_SpAMM_RNK2, &
        SpAMM_New_SpAMM_C, &
        SpAMM_New_QuTree

  END INTERFACE SpAMM_New

  INTERFACE SpAMM_Set

    MODULE PROCEDURE &
        SpAMM_Set_SpAMM_RNK2

  END INTERFACE SpAMM_Set

  INTERFACE SpAMM_Get

    MODULE PROCEDURE &
        SpAMM_Get_SpAMM_RNK2, &
        SpAMM_Get_SpAMM_C

  END INTERFACE SpAMM_Get

  INTERFACE SpAMM_SetEq

    MODULE PROCEDURE &
        SpAMM_SetEq_SpAMM_RNK2_to_Dense, &
        SpAMM_SetEq_SpAMM_C_to_Dense, &
        SpAMM_SetEq_SpAMM_C_to_SpAMM_C, &
        SpAMM_SetEq_Dense_to_SpAMM_C

  END INTERFACE SpAMM_SetEq

  INTERFACE

    !> Interface for spamm_check().
    SUBROUTINE spamm_check (A) &
        BIND(C, NAME = "spamm_check_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr), INTENT(IN) :: A
    END SUBROUTINE spamm_check

    !> Interface for spamm_get_norm().
    FUNCTION spamm_get_norm (A) &
        BIND(C, NAME = "spamm_get_norm_interface")
      USE, INTRINSIC :: iso_C_binding
      USE SpAMMPACK_TYPES
      REAL(NORM_TYPE) :: spamm_get_norm
      TYPE(c_ptr), INTENT(IN) :: A
    END FUNCTION spamm_get_norm

    !> Interface for spamm_get_number_dimensions().
    FUNCTION spamm_get_number_dimensions (A) &
        BIND(C, NAME = "spamm_get_number_dimensions_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int) :: spamm_get_number_dimensions
      TYPE(c_ptr), INTENT(IN) :: A
    END FUNCTION spamm_get_number_dimensions

    !> Interface for spamm_get_N().
    SUBROUTINE spamm_get_N (M, N, A) &
        BIND(C, NAME = "spamm_get_N_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int), INTENT(INOUT) :: M, N
      TYPE(c_ptr), INTENT(IN) :: A
    END SUBROUTINE spamm_get_N

    !> Interface for spamm_convert_dense_to_spamm().
    SUBROUTINE spamm_convert_dense_to_spamm (ndim, N, chunk_tier, use_linear_tree, A_dense, A) &
        BIND(C, NAME = "spamm_convert_dense_to_spamm_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int), INTENT(IN)                  :: ndim
      INTEGER(c_int), DIMENSION(ndim), INTENT(IN) :: N
      INTEGER(c_int), INTENT(IN)                  :: chunk_tier
      INTEGER(c_int), INTENT(IN)                  :: use_linear_tree
      REAL(c_float), DIMENSION(*), INTENT(IN)     :: A_dense
      TYPE(c_ptr), INTENT(INOUT)                  :: A
    END subroutine spamm_convert_dense_to_spamm

    !> Interface for spamm_convert_spamm_to_dense().
    SUBROUTINE spamm_convert_spamm_to_dense (A, ADense) &
        BIND(C, name = "spamm_convert_spamm_to_dense_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr), INTENT(IN)                 :: A
      REAL(c_float), DIMENSION(*), INTENT(IN) :: ADense
    END SUBROUTINE spamm_convert_spamm_to_dense

    !> Interface for spamm_get().
    SUBROUTINE spamm_get_element (i, A, Aij) &
        BIND(C, name = "spamm_get_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int), DIMENSION(*) :: i
      TYPE(c_ptr), INTENT(IN)      :: A
      REAL(c_float), INTENT(OUT)  :: Aij
    END SUBROUTINE spamm_get_element

    SUBROUTINE spamm_new_interface (ndim, N, chunkTier, useLinearTree, A) &
        BIND(C, name = "spamm_new_interface")
      USE, INTRINSIC :: iso_C_binding
      INTEGER(c_int), INTENT(IN)               :: ndim
      INTEGER(c_int), DIMENSION(2), INTENT(IN) :: N
      INTEGER(c_int), INTENT(IN)               :: chunkTier
      INTEGER(c_int), INTENT(IN)               :: useLinearTree
      TYPE(c_ptr), INTENT(INOUT)               :: A
    END SUBROUTINE spamm_new_interface

    SUBROUTINE spamm_copy (A, beta, B, flop, mop) &
        BIND(C, name = "spamm_copy_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr), INTENT(INOUT)    :: A
      REAL(c_float), INTENT(IN)     :: beta
      TYPE(c_ptr), INTENT(IN)       :: B
      REAL(c_double), INTENT(INOUT) :: flop
      REAL(c_double), INTENT(INOUT) :: mop
    END SUBROUTINE spamm_copy

    SUBROUTINE spamm_print (A) &
        BIND(C, name = "spamm_print_interface")
      USE, INTRINSIC :: iso_C_binding
      TYPE(c_ptr), INTENT(IN) :: A
    END SUBROUTINE spamm_print

  END INTERFACE

CONTAINS

  SUBROUTINE SpAMM_New_QuTree (qA)

    TYPE(QuTree), INTENT(INOUT) :: qA

  END SUBROUTINE SpAMM_New_QuTree

  INTEGER FUNCTION get_tree_depth (numberDimensions, N, useLinearTree) RESULT(depth)

    INTEGER, INTENT(IN)                              :: numberDimensions
    INTEGER, DIMENSION(numberDimensions), INTENT(IN) :: N
    LOGICAL, INTENT(IN)                              :: useLinearTree

    INTEGER :: NTemp, dim
    REAL*4 :: x, xN

    x = 0E0;
    DO dim = 1, numberDimensions
      IF(useLinearTree) THEN
        IF(N(dim) < SPAMM_N_KERNEL) THEN
          NTemp = SPAMM_N_KERNEL
        ELSE
          NTemp = N(dim)
        ENDIF
      ELSE
        NTemp = N(dim)
      ENDIF

      xN = log(REAL(NTemp))/log(2E0)

      IF(xN > x) THEN
        x = xN
      ENDIF
    ENDDO

    depth = ceiling(x)

    ! Double check the depth.
    IF(depth >= 1) THEN
      DO dim = 1, numberDimensions
        IF(2**(depth-1) < N(dim)) THEN
          depth = depth+1
          EXIT
        ENDIF
      ENDDO
      depth = depth-1
    ENDIF

  END FUNCTION get_tree_depth

  SUBROUTINE SpAMM_New_SpAMM_RNK2 (N, chunkTier, useLinearTree, A)

    INTEGER, DIMENSION(2), INTENT(IN) :: N
    INTEGER, INTENT(IN)               :: chunkTier
    LOGICAL, INTENT(IN)               :: useLinearTree
    TYPE(SpAMM_RNK2), INTENT(INOUT)   :: A

    A%N = N
    A%depth = get_tree_depth(2, N, useLinearTree)
    A%NPadded = 2**A%depth
    IF(chunkTier > A%depth) THEN
      A%chunkTier = A%depth
    ELSE
      A%chunkTier = chunkTier
      A%depth = chunkTier
    ENDIF

  END SUBROUTINE SpAMM_New_SpAMM_RNK2

  SUBROUTINE SpAMM_New_SpAMM_C (N, chunkTier, useLinearTree, A)

    INTEGER, DIMENSION(2), INTENT(IN) :: N
    INTEGER, INTENT(IN)               :: chunkTier
    LOGICAL, INTENT(IN)               :: useLinearTree
    TYPE(c_ptr), INTENT(INOUT)        :: A

    INTEGER(c_int) :: use_linear_tree

    IF(useLinearTree) THEN
      use_linear_tree = 1
    ELSE
      use_linear_tree = 0
    ENDIF

    CALL spamm_new_interface(2, N, chunkTier, use_linear_tree, A)

  END SUBROUTINE SpAMM_New_SpAMM_C

  SUBROUTINE SpAMM_Set_SpAMM_RNK2 (i, j, Aij, A)

    INTEGER, INTENT(IN)             :: i, j
    REAL*4, INTENT(IN)              :: Aij
    TYPE(SpAMM_RNK2), INTENT(INOUT) :: A

  END SUBROUTINE SpAMM_Set_SpAMM_RNK2

  REAL*4 FUNCTION SpAMM_Get_SpAMM_RNK2 (i, j, A) RESULT(Aij)

    INTEGER, INTENT(IN)          :: i, j
    TYPE(SpAMM_RNK2), INTENT(IN) :: A

    Aij = 0

  END FUNCTION SpAMM_Get_SpAMM_RNK2

  REAL*4 FUNCTION SpAMM_Get_SpAMM_C (i, j, A) RESULT(Aij)

    INTEGER, INTENT(IN) :: i, j
    TYPE(c_ptr), INTENT(IN) :: A

    CALL spamm_get_element((/ i-1, j-1 /), A, Aij)

  END FUNCTION SpAMM_Get_SpAMM_C

  SUBROUTINE SpAMM_SetEq_SpAMM_RNK2_to_Dense (A, ADense)

    TYPE(SpAMM_RNK2), INTENT(INOUT)    :: A
    REAL*4, DIMENSION(:,:), INTENT(IN) :: ADense

    INTEGER :: i, j

    IF(A%N(1) /= SIZE(ADense, 1)) THEN
      STOP
    ENDIF

    IF(A%N(2) /= SIZE(ADense, 2)) THEN
      STOP
    ENDIF

    DO i = 1, A%N(1)
      DO j = 1, A%N(2)
        CALL SpAMM_Set(i, j, ADense(i,j), A)
      ENDDO
    ENDDO

  END SUBROUTINE SpAMM_SetEq_SpAMM_RNK2_to_Dense

  !> Convert a dense matrix to a SpAMM matrix.
  !!
  !! @param A The SpAMM matrix.
  !! @param ADense The dense matrix.
  !! @param chunkTier The chunk tier to use for the SpAMM matrix.
  !! @param useLinearTree Whether to use a linear tree at the chunk tier or not.
  SUBROUTINE SpAMM_SetEq_SpAMM_C_to_Dense (A, ADense, chunkTier, useLinearTree)

    TYPE(c_ptr), INTENT(INOUT)         :: A
    REAL*8, DIMENSION(:,:), INTENT(IN) :: ADense
    INTEGER, INTENT(IN)                :: chunkTier
    LOGICAL, INTENT(IN)                :: useLinearTree

    REAL*4, DIMENSION(SIZE(ADense, 1), SIZE(ADense, 2)) :: ADenseSingle
    INTEGER :: i, j
    INTEGER :: use_linear_tree

    IF(useLinearTree) THEN
      use_linear_tree = 1
    ELSE
      use_linear_tree = 0
    ENDIF

    DO i = 1, size(ADense, 1)
      DO j = 1, size(ADense, 2)
        ADenseSingle(i,j) = ADense(i,j)
      ENDDO
    ENDDO

    CALL spamm_convert_dense_to_spamm(2, (/ size(ADense, 1), size(ADense, 2) /), &
      chunkTier, use_linear_tree, ADenseSingle, A)

  END SUBROUTINE SpAMM_SetEq_SpAMM_C_to_Dense

  SUBROUTINE SpAMM_SetEq_SpAMM_C_to_SpAMM_C (A, B, flop_O, mop_O)

    TYPE(c_ptr), INTENT(INOUT)      :: A
    TYPE(c_ptr), INTENT(IN)         :: B
    REAL*8, INTENT(INOUT), OPTIONAL :: flop_O
    REAL*8, INTENT(INOUT), OPTIONAL :: mop_O

    REAL(c_double) :: flop
    REAL(c_double) :: mop

    flop = 0.0D0
    mop = 0.0D0

    CALL spamm_copy(A, 1.0, B, flop, mop)

    IF(PRESENT(flop_O)) THEN
      flop_O = flop_O+flop
    ENDIF
    IF(PRESENT(mop_O)) THEN
      mop_O = mop_O+mop
    ENDIF

  END SUBROUTINE SpAMM_SetEq_SpAMM_C_to_SpAMM_C

  !> Convert a SpAMM matrix to a dense matrix.
  !!
  !! @param ADense The dense matrix.
  !! @param A The SpAMM matrix.
  SUBROUTINE SpAMM_SetEq_Dense_to_SpAMM_C (ADense, A)

    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: ADense
    TYPE(c_ptr), INTENT(IN)               :: A

    REAL*4 :: temp

    INTEGER :: i, j
    INTEGER :: M, N

    ! Get matrix size.
    CALL spamm_get_N(M, N, A)

    IF(size(ADense, 1) /= M .OR. size(ADense, 2) /= N) THEN
      WRITE(*, *) "dimension mismatch, ADense = ", size(ADense, 1), size(ADense, 2)
      WRITE(*, *) "dimension mismatch, SpAMM = ", M, N
      STOP
    ENDIF

    DO i = 1, M
      DO j = 1, N
        temp = SpAMM_Get(i, j, A)
        ADense(i, j) = temp

        ! This is an ugly hack to fix an apparent compiler bug. Without those
        ! lines the compiler translates the 2 lines above as:
        !
        ! 356	        temp = SpAMM_Get(i, j, A)
        !    0x000000000000046b <+1131>:	mov    -0x4d0(%rbp),%rdx
        !    0x0000000000000472 <+1138>:	lea    -0x54(%rbp),%rcx
        !    0x0000000000000476 <+1142>:	lea    -0x50(%rbp),%rax
        !    0x000000000000047a <+1146>:	mov    %rcx,%rsi
        !    0x000000000000047d <+1149>:	mov    %rax,%rdi
        !    0x0000000000000480 <+1152>:	callq  0x485 <__spammpack_management_MOD_spamm_seteq_dense_to_spamm_c+1157>
        !    0x0000000000000485 <+1157>:	movss  %xmm0,-0x4c(%rbp)
        !
        ! 357	        ADense(i, j) = temp
        !    0x000000000000048a <+1162>:	mov    -0x50(%rbp),%eax
        !    0x000000000000048d <+1165>:	cltq
        !    0x000000000000048f <+1167>:	mov    %rax,%rdx
        !    0x0000000000000492 <+1170>:	imul   %rbx,%rdx
        !    0x0000000000000496 <+1174>:	mov    -0x54(%rbp),%eax
        !    0x0000000000000499 <+1177>:	cltq
        !    0x000000000000049b <+1179>:	imul   %r12,%rax
        !    0x000000000000049f <+1183>:	add    %rdx,%rax
        !    0x00000000000004a2 <+1186>:	lea    (%rax,%r15,1),%rdx
        !    0x00000000000004a6 <+1190>:	movss  -0x4c(%rbp),%xmm0
        !    0x00000000000004ab <+1195>:	cvtps2pd %xmm0,%xmm0
        !    0x00000000000004ae <+1198>:	mov    -0x38(%rbp),%rax
        !    0x00000000000004b2 <+1202>:	movsd  %xmm0,(%rax,%rdx,8)
        !
        ! With the two lines below, the compiler inserts an additional
        ! statement:
        !
        !    0x00000000000004af <+1199>:	unpcklps %xmm0,%xmm0
        !
        ! When spamm_get() returns a zero, the high quadword in the double after
        ! cvtps2pd contains junk and ADense ends up with very large numbers.
        IF(abs(temp) > 1.0e4) THEN
          WRITE(*, *) "large matrix element: ", i, j, temp
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE SpAMM_SetEq_Dense_to_SpAMM_C

END MODULE SpAMMPACK_MANAGEMENT
