!> The SpAMM data structures, STDEC and SALDS
module spamm_structures

  use spamm_parameters

  implicit none

  ! SpAMM tree decorations  ________________ STDECs _________________

  !> garnishments of the tree_1d ...
  type :: SpAMM_decoration_1d
     !> initialization status
     logical :: is_initialized = .false.
     !> leaf node flag
     logical :: leaf = .false.
     !> tree sub-vector width: M_pad/2**depth
     integer :: width = -1
     !> Integer dimension of the native (non-padded) vector
     integer :: ndimn = -1
     !> Axis-aligned bounding box for the [i] index space
     integer,  dimension(0:1) :: BndBx
     !> Square of the F-norm.
     real(SPAMM_KIND) :: Norm2 = -1
     !> Float Ops needed accumulated to this level
     real(kind(0d0)) :: FlOps = -1
     !> The number of non-zero elements to this level
     real(kind(0d0)) :: Non0s = -1
  end type SpAMM_decoration_1d

  !> Garnishments of tree_2d ...
  type :: SpAMM_decoration_2d
     !> initialization status
     logical :: is_initialized = .false.
     !> leaf node flag
     logical :: leaf = .false.
     !> tree sub-matrix width: MN_pad/2**depth
     integer, dimension(1:2) :: width = -1
     !> Integer dimension of the native (non-padded) matrix
     integer, dimension(1:2) :: ndimn = -1
     !> Axis-aligned bounding box for the [i]-[j] index space
     integer, dimension(0:1, 1:2) :: bndbx
     !> Square of the F-norm.
     real(SPAMM_KIND) :: norm2 = -1
     !> Float Ops needed accumulated to this level
     real(kind(0d0)) :: flops = -1
     !> The number of non-zero elements to this level
     real(kind(0d0)) :: non0s = -1
  end type SpAMM_decoration_2d

  ! SpAMM algebraic data structures _______________ SALGDSs ____________________

  !> The tree_1d (vector) type:
  type :: SpAMM_tree_1d
     type(SpAMM_decoration_1d)             :: frill
     type(SpAMM_tree_1d),      pointer     :: child_0 => null()
     type(SpAMM_tree_1d),      pointer     :: child_1 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:)
  end type SpAMM_tree_1d

  !> The tree_2d matrix structures:
  !! symmetric (SPD/Hermetian) ...
  type :: SpAMM_tree_2d_symm
     type(SpAMM_decoration_2d)             :: frill
     type(SpAMM_tree_2d_symm), pointer     :: child_00 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_01 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_10 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_11 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_2d_symm

  !> full ...
  type :: SpAMM_tree_2d_full
     type(SpAMM_decoration_2d)             :: frill
     type(SpAMM_tree_2d_full), pointer     :: child_00 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_01 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_10 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_11 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_2d_full

end module spamm_structures
