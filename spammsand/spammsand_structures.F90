!----------------------------------------------------------------------------------
! Nested SpAMM solvers (SpAMM sandwitches) for matrix functions
MODULE SpAMMsand_structures

  USE  spammpack

  implicit none

  ! SpAMM sandwich for matrix functions: f(|a>) = |x_0>.|x_1> ... |x_m>
  type :: spammsand_tree_2d_slices

     ! the SpAMM squared threshold for each algebra
     real(SpAMM_KIND)                        :: tau_0
     real(SpAMM_KIND)                        :: tau_S
     ! the first approximation its nested residuals ...
     type(SpAMM_tree_2d_symm),       pointer :: mtx
     type(spammsand_tree_2d_slices), pointer :: nxt

  end type spammsand_tree_2d_slices

  ! --
contains

  ! --
end module spammsand_structures ! ... and we're out ...
