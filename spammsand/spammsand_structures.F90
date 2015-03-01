!----------------------------------------------------------------------------------
! Nested SpAMM solvers (SpAMM sandwitches) for matrix functions
MODULE SpAMMsand_structures

  USE  spammpack

  implicit none

  integer, parameter :: slices=4

  ! SpAMM sandwich 4 matrix functions: f(|a>) = |x_0>.|x_1> ... |x_m>
  type :: spammsand_tree_2d_symm

     ! the SpAMM squared threshold for each algebra
     integer,                  dimension(:), pointer  :: tau2
     ! the factor and its nested residuals, a sandwich ...
     type(SpAMM_tree_2d_symm), dimension(:), pointer  :: trix

  end type spammsand_tree_2d_symm

  ! --
contains

  ! --
end module spammsand_structures ! ... and we're out ...
