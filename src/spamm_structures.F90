!> The SpAMM data structures, STDEC and SALDS
module spamm_structures

  use spamm_parameters

  implicit none

  !> Garnishments of the tree_1d ...
  type :: spamm_decoration_1d
     !> initialization status
     logical :: needs_initialization = .false.
     !> leaf node flag
     logical :: leaf = .false.
     !> tree sub-vector width: M_pad/2**depth
     integer :: width = -1
     !> Integer dimension of the native (non-padded) vector
     integer :: ndimn = -1
     !> Axis-aligned bounding box for the [i] index space in the
     !> unpadded dimensions.
     integer,  dimension(0:1) :: bndbx
     !> Square of the F-norm.
     real(SPAMM_KIND) :: norm2 = -1
     !> Float Ops needed accumulated to this level
     double precision :: flops = -1
     !> The number of non-zero elements to this level
     double precision :: non0s = -1
  end type spamm_decoration_1d

  !> Garnishments of tree_2d ...
  type :: spamm_decoration_2d
     !> initialization status
     logical :: needs_initialization = .false.
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
  end type spamm_decoration_2d

  !> The tree_1d (vector) type.
  !! @ingroup types_group
  type :: spamm_tree_1d
     !> Tree decoration.
     type(spamm_decoration_1d) :: frill
     !> Pointer to quadrant.
     type(spamm_tree_1d), pointer :: child_0 => null()
     !> Pointer to quadrant.
     type(spamm_tree_1d), pointer :: child_1 => null()
     !> Matrix data at leaf nodes.
     real(SPAMM_KIND), allocatable :: chunk(:)
  end type spamm_tree_1d

  !> The tree_2d matrix structures.
  !! symmetric (SPD/Hermetian)
  !! @ingroup types_group
  type :: spamm_tree_2d_symm
     !> Tree decoration.
     type(spamm_decoration_2d) :: frill
     !> Pointer to quadrant.
     type(spamm_tree_2d_symm), pointer :: child_00 => null()
     !> Pointer to quadrant.
     type(spamm_tree_2d_symm), pointer :: child_01 => null()
     !> Pointer to quadrant.
     type(spamm_tree_2d_symm), pointer :: child_10 => null()
     !> Pointer to quadrant.
     type(spamm_tree_2d_symm), pointer :: child_11 => null()
     !> Matrix data at leaf nodes.
     real(SPAMM_KIND), allocatable :: chunk(:, :)
  end type spamm_tree_2d_symm

  !> full ...
  !! @ingroup types_group
  type :: spamm_tree_2d_full
     !> Tree decoration.
     type(spamm_decoration_2d) :: frill
     !> Pointer to quadrant.
     type(spamm_tree_2d_full), pointer :: child_00 => null()
     !> Pointer to quadrant.
     type(spamm_tree_2d_full), pointer :: child_01 => null()
     !> Pointer to quadrant.
     type(spamm_tree_2d_full), pointer :: child_10 => null()
     !> Pointer to quadrant.
     type(spamm_tree_2d_full), pointer :: child_11 => null()
     !> Matrix data at leaf nodes.
     real(SPAMM_KIND), allocatable :: chunk(:, :)
  end type spamm_tree_2d_full

end module spamm_structures
