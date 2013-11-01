MODULE SpAMMPACK_TYPES

  USE, INTRINSIC :: iso_C_binding

#ifdef _OPENMP
  USE omp_lib
#endif

  IMPLICIT NONE

  !> Type for norms and tolerances.
#if SPAMM_NORM_TYPE == float
  INTEGER, PARAMETER :: NORM_TYPE = c_float
#elif SPAMM_NORM_TYPE == double
  INTEGER, PARAMETER :: NORM_TYPE = c_double
#else
#error "unknown type"
#endif

  !> Rank 2 SpAMM matrix.
  TYPE SpAMM_RNK2

    !> Number of rows and columns.
    INTEGER, DIMENSION(2) :: N

    !> Padded matrix size.
    INTEGER :: NPadded

    !> The total depth of the tree.
    INTEGER :: depth = -1

    !> The tier of the contiguous matrix chunks. At this tier the submatrix is
    !> stored in SpAMM chunks. If use_linear_tree, then the spamm_hashed_* code
    !> takes over and processes the chunks, if not, then the chunk is
    !> interpreted has containing a submatrix of dimension N_contiguous x
    !> N_contiguous. */
    INTEGER :: chunkTier

    !> Use a linear submatrix representation or a hierarchical one. A value of
    !> 0 implies hierarchical, while anything else (but we use 1) implies
    !> linear. */
    LOGICAL :: useLinearTree

    !> The tree.
    TYPE(QuTree), POINTER :: root => NULL()

  END TYPE SpAMM_RNK2

  !> A pointer to a Quadtree object.
  TYPE QuTreePointer

    !> The pointer.
    TYPE(QuTree), POINTER :: node => NULL()

  END TYPE QuTreePointer

  !> Quaternary tree data structure.
  TYPE QuTree

    !> The Frobenious norm.
    REAL*4 :: norm

    !> The squared Frobenious norm.
    REAL*4 :: norm2

    !> The pointers to the subtrees.
    TYPE(QuTreePointer), DIMENSION(2,2) :: child

    !> The matrix data.
    TYPE(c_ptr) :: chunk

#ifdef _OPENMP
    !> Block lock
    INTEGER(KIND = OMP_LOCK_KIND) :: lock
#endif

  END TYPE QuTree

END MODULE SpAMMPACK_TYPES
