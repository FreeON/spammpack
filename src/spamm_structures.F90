!----------------------------------------------------------------------------------
! The SpAMM data structures, STDEC and SALDS
!
MODULE spamm_structures

  USE  spamm_parameters

  implicit none

  ! SpAMM tree decorations  ________________ STDECs _________________

  ! garnishments of the tree_1d ...
  type :: SpAMM_decoration_1d
     ! initialization status
     logical                               :: Init
     ! leaf node flag
     logical                               :: Leaf
     !> tree sub-vector width: M_pad/2**depth
     integer                               :: Width
     !> Integer dimension of the native (non-padded) vector
     integer                               :: NDimn
     !> Axis-aligned bounding box for the [i] index space
     integer,  dimension(0:1)              :: BndBx
     !> Square of the F-norm.
     real(SPAMM_KIND)                      :: Norm2 = -1
     !> Float Ops needed accumulated to this level
     real(kind(0d0))                       :: FlOps = -1
     !> The number of non-zero elements to this level
     real(kind(0d0))                       :: Non0s = -1
  end type SpAMM_decoration_1d

  ! garnishments of tree_2d ...
  type :: SpAMM_decoration_2d
     ! initialization status
     logical                               :: Init
     ! leaf node flag
     logical                               :: Leaf
     !> tree sub-matrix width: MN_pad/2**depth
     integer,           dimension(1:2)     :: Width
     !> Integer dimension of the native (non-padded) matrix
     integer,           dimension(1:2)     :: NDimn
     !> Axis-aligned bounding box for the [i]-[j] index space
     integer,  dimension(0:1,1:2)          :: BndBx
          !> Square of the F-norm.
     real(SPAMM_KIND)                      :: Norm2 = -1
     !> Float Ops needed accumulated to this level
     real(kind(0d0))                       :: FlOps = -1
     !> The number of non-zero elements to this level
     real(kind(0d0))                       :: Non0s = -1
  end type SpAMM_decoration_2d

  ! SpAMM algebraic data structures _______________ SALGDSs ____________________

  ! The tree_1d (vector) type:
  type :: SpAMM_tree_1d
     type(SpAMM_decoration_1d)             :: frill
     type(SpAMM_tree_1d),      pointer     :: child_0 => null()
     type(SpAMM_tree_1d),      pointer     :: child_1 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:)
  end type SpAMM_tree_1d

  ! The tree_2d matrix structures:
  ! symmetric (SPD/Hermetian) ...
  type :: SpAMM_tree_2d_symm
     type(SpAMM_decoration_2d)             :: frill
     type(SpAMM_tree_2d_symm), pointer     :: child_00 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_01 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_10 => null()
     type(SpAMM_tree_2d_symm), pointer     :: child_11 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_2d_symm

  ! full ...
  type :: SpAMM_tree_2d_full
     type(SpAMM_decoration_2d)             :: frill
     type(SpAMM_tree_2d_full), pointer     :: child_00 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_01 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_10 => null()
     type(SpAMM_tree_2d_full), pointer     :: child_11 => null()
     real(SPAMM_KIND),     allocatable     :: chunk(:, :)
  end type SpAMM_tree_2d_full

  ! --
contains

  ! --
end module spamm_structures ! ... and we're out ...
